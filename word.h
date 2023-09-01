#ifndef WORD_H
#define WORD_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "mt64.h"

#ifndef UINT64_MAX
#error No 64-bit unsigned integer type!
#endif

typedef uint64_t word_t;
#define DECWLEN (20)
#define PRIw PRIu64

#define WBITS  ((int)64)
#define WBITS1 ((int)63)

#define WONE   ((word_t)1)
#define WZERO  ((word_t)0)
#define WONES  (~WZERO)
#define WHIONE (WONE<<WBITS1)

#define ZWORD(w)  ((w) == WZERO)
#define NZWORD(w) ((w) != WZERO)

#define POW2(n) (WONE << (n))

#define WHMAX (WONE << 32)

// Return WONE or WZERO (true or false) depending whether bit at position p in w is set; w should be of type word_t
#define BITON(w,p) (WONE&((w)>>(p)))

// Get bit in word w at position p (i.e. mask other bits to 0); w should be of type word_t
#define GETBIT(w,p) ((WONE<<(p))&(w))

// Set bit in word w at position p; w should be of type word_t
#define SETBIT(w,p) ((w) |= (WONE<<(p)))

// Set bit in word w at position p to bit b; w should be of type word_t, and b must be either WONE or WZERO
// Note that if p is an expression, it will be evaluated twice in undefined order - careful!
#define PUTBIT(w,p,b) ((w) = (((w)&(~(WONE<<(p))))|((b)<<(p))))

// Flip bit in word w at position p; w should be of type word_t
#define FLIPBIT(w,p) ((w) ^= (WONE<<(p)))

// Usual caveats!
#define SWAP(type,a,b) {type const a##tmp = a; a = b; b = a##tmp;}

// other handy stuff
#define newline putchar('\n')
#define fnewline(fstream) fputc('\n',fstream)

/*********************************************************************/
/*                      single-word                                  */
/*********************************************************************/

static inline word_t wd_random(mt_t* const prng)
{
	return mt_uint(prng);
}

static inline word_t wd_randomb(const double p, mt_t* const prng)
{
	word_t w = WZERO;
	for (int i=0; i<WBITS; ++i) if (mt_rand(prng) < p) SETBIT(w,i);
	return w;
}

static inline int wd_cointoss(mt_t* const prng)
{
	return mt_uint(prng) < WHMAX;
}

static inline word_t wd_rotl(const word_t w, const int b)
{
	return (w<<b)|(w>>(WBITS-b));
}

static inline word_t wd_rotr(const word_t w, const int b)
{
	return (w>>b)|(w<<(WBITS-b));
}

static inline size_t wd_filter(const int n, const word_t w, const int B, const word_t* const f)
{
	// CAUTION: need 2*n <= WBITS, and assumes WBITS-n hi-bits of w are cleared!!!
	const word_t BMASK = WONES>>(WBITS-B); // mask to clear bits above 1st B
	word_t w2 = (w<<n)|w; // double-up word
	word_t wnew = WZERO;
	for (int i=0; i<n; ++i, w2>>=1) PUTBIT(wnew,i,f[w2&BMASK]);
	return wnew;
}

static inline word_t wd_reverse(word_t w)
// See http://graphics.stanford.edu/~seander/bithacks.html for even funkier versions.
{
	int s = (int)WBITS1; // extra shift needed at end
	for (word_t u = w >> 1; u; u >>= 1,--s) {
		w <<= 1;
		w |= u&1;
	}
	w <<= s; // shift when u's highest bits are zero
	return w;
}

static inline void wd_noisify(word_t* const w, const double p, mt_t* const prng)
{
	for (int i=0; i<WBITS; ++i) if (mt_rand(prng) < p) FLIPBIT(*w,i);
}

static inline int wd_nsetbits(word_t w) // Kernighan!
{
// See http://graphics.stanford.edu/~seander/bithacks.html for even funkier versions.
	int c; // c accumulates the total bits set in w
	for (c = 0;w;++c) w &= w-1; // clear the least significant bit set
	return c;
}

void wd_print    (const word_t w);
void wd_fprint   (const word_t w, FILE* const fstream);
void wd_prints   (const word_t w);
void wd_fprints  (const word_t w, FILE* const fstream);
void wd_printc   (const word_t w, const int W);
void wd_fprintc  (const word_t w, const int W, FILE* const fstream);
void wd_print_lo (const word_t w, const int b);
void wd_fprint_lo(const word_t w, const int b, FILE* const fstream);

void wd_dft(const word_t w, double* const wdftre, double* const wdftim, const double* const costab, const double* const sintab);

void wd_autocov(const word_t w, double* const wac);

/*********************************************************************/
/*                      multi-word                                   */
/*********************************************************************/

word_t* mw_alloc(const size_t n);

word_t* mw_alloc_copy(const size_t n, const word_t* const wsrc);

static inline void mw_copy(const size_t n, word_t* const wdest, const word_t* const wsrc)
{
	// NOTE: wdest and wsrc must not overlap!!!
	memcpy(wdest,wsrc,n*sizeof(word_t));
}

static inline void mw_zero(const size_t n, word_t* const w)
{
	memset(w,0,n*sizeof(word_t));
}

static inline void mw_randomise(const size_t n, word_t* const w, mt_t* const prng)
{
	for (word_t* pw=w;pw<w+n;++pw) *pw = wd_random(prng);
}

static inline void mw_randomiseb(const size_t n, word_t* const w, const double p, mt_t* const prng)
{
	for (word_t* pw=w;pw<w+n;++pw) *pw = wd_randomb(p,prng);
}

static inline int mw_equal(const size_t n, const word_t* const w1, const word_t* const w2)
{
	return memcmp(w1,w2,n*sizeof(word_t)) == 0;
//	for (size_t k=0;k<n;++k) if (w2[k] != w1[k]) return 0;
//	return 1;
}

static inline void mw_rotl(const size_t n, word_t* const wrot, const word_t* const w, const size_t nbits)
{
	const int b  = nbits%WBITS;
	const int Wb = WBITS-b;
	const size_t m = (nbits/WBITS)%n;
	if (b == 0) {
		// could use memmove() here, but maybe faster without unless n is quite large
		size_t k = 0;
		for (;k<n-m;++k) wrot[k+m  ] = w[k];
		for (;k<n;  ++k) wrot[k+m-n] = w[k];
		return;
	}
	wrot[m] = (w[0]<<b)|(w[n-1]>>Wb);
	size_t k = 1;
	for (;k<n-m;++k) wrot[k+m  ] = (w[k]<<b)|(w[k-1]>>Wb);
	for (;k<n;  ++k) wrot[k+m-n] = (w[k]<<b)|(w[k-1]>>Wb);
}

static inline void mw_rotr(const size_t n, word_t* const wrot, const word_t* const w, const size_t nbits)
{
	mw_rotl(n,wrot,w,n*WBITS-nbits); // go the long way round :-)
}

static inline void mw_reverse(const size_t n, word_t* const wrev, const word_t* const w)
{
	size_t n1 = n-1;
	for (size_t k=0; k<n; ++k) wrev[n1-k] = wd_reverse(w[k]);
}

static inline void wm_noisify(const size_t n, word_t* const w, const double p, mt_t* const prng)
{
	for (word_t* pw=w; pw<w+n; ++pw) wd_noisify(pw,p,prng);
}

static inline int mw_equiv(const size_t n, const word_t* const w1, const word_t* const w2)
{
	word_t w2rot[n];
	for (size_t b=0;b<n*WBITS;++b) {
		mw_rotl(n,w2rot,w2,b);
		if (mw_equal(n,w1,w2rot)) return (int)b;
	}
	return -1;
}

static inline int mw_nsetbits(const size_t n, const word_t* const w)
{
	int b = 0;
	for (const word_t* u=w;u<w+n;++u)  b += wd_nsetbits(*u);
	return b;
}

static inline int mw_iszero(const size_t n, const word_t* const w)
{
	for (const word_t* u=w;u<w+n;++u) if (*u) return 0;
	return 1;
}

static inline void mw_filter(const size_t n, word_t* const wnew, const word_t* const w, const int B, const word_t* const f)
{
	const int WB = WBITS-B;
	const word_t BMASK = WONES>>WB; // mask to clear bits above 1st B
	for (size_t k=0;k<n;++k) {
		size_t wk = w[k];
		const word_t wk1 = k < n-1 ? w[k+1] : w[0]; // next word : wrap to lo-word on last word
		int i = 0;
		word_t wnewk;
		for (;i<WB;   ++i,wk>>=1) PUTBIT(wnewk,i,f[wk&BMASK]);
		wk |= (wk1<<B); // splice in next word
		for (;i<WBITS;++i,wk>>=1) PUTBIT(wnewk,i,f[wk&BMASK]);
		wnew[k] = wnewk;
	}
}

void mw_run(const size_t I, const size_t n, word_t* const w, const int B, const word_t* const f);

static inline word_t mw_get_part(const size_t n, const word_t* const w, const int B, const int b)
{
	const int WB = WBITS-B;
	const int WE = WBITS+WB;
	const size_t k = (size_t)(b/WBITS);
	const int i = b%WBITS;
	return (i<=WB?(w[k]<<(WB-i))>>WB:(w[k]>>i)|((w[k<n-1?k+1:0]<<(WE-i))>>WB));
}

static inline void mw_parts(const size_t n, word_t* const wnew, const word_t* const w, const int B, const word_t part)
{
	const int WB = WBITS-B;
	const int WE = WBITS+WB;
	size_t k = 0;
	for (;k<n-1;++k) {
		int i = 0;
		for (;i<=WB;  ++i) PUTBIT(wnew[k],i,(word_t)(((w[k]<<(WB-i))>>WB) == part));
		for (;i<WBITS;++i) PUTBIT(wnew[k],i,(word_t)(((w[k]>>i)|((w[k+1]<<(WE-i))>>WB)) == part));
	}
	{ // k = n-1 wrap lo-word
		int i = 0;
		for (;i<=WB;  ++i) PUTBIT(wnew[k],i,(word_t)(((w[k]<<(WB-i))>>WB) == part));
		for (;i<WBITS;++i) PUTBIT(wnew[k],i,(word_t)(((w[k]>>i)|((w[0]<<(WE-i))>>WB)) == part));
	}
}

typedef struct {word_t w; double f;} partf_t;
static inline int partf_comp(const void* const elem1, const void* const elem2)
{
    const double f1 = ((partf_t*)elem1)->f;
    const double f2 = ((partf_t*)elem2)->f;
    return (f2 > f1 ? 1 : f2 < f1 ? -1 : 0);
}

void mw_fprint     (const size_t n, const word_t* const w, FILE* const fstream);
void mw_print      (const size_t n, const word_t* const w);
void mw_fprints    (const size_t n, const word_t* const w, FILE* const fstream);
void mw_prints     (const size_t n, const word_t* const w);
void mw_fprint_bin (const size_t n, const word_t* const w, FILE* const fstream);
void mw_print_bin  (const size_t n, const word_t* const w);

void mw_dft(const size_t n, const word_t* const w, double* const dftre, double* const dftim, double* const dsp, const double* const costab);

void mw_autocov(const size_t n, const word_t* const w, double* const ac);

/*********************************************************************/

// printf binary macros

#define PBP08 "%c%c%c%c%c%c%c%c"
#define PBP16 PBP08 PBP08
#define PBP32 PBP16 PBP16
#define PBP64 PBP32 PBP32
#define PBP PBP64

#define PBC08 PBP08
#define PBC16 PBC08","PBC08
#define PBC32 PBC16","PBC16
#define PBC64 PBC32","PBC32
#define PBC PBC64

#define PBS08 " %c %c %c %c %c %c %c %c"
#define PBS16 PBS08 PBS08
#define PBS32 PBS16 PBS16
#define PBS64 PBS32 PBS32
#define PBS PBS64

#define PBI08(i) \
    (((i) & 0x80ll) ? '1' : '0'), \
    (((i) & 0x40ll) ? '1' : '0'), \
    (((i) & 0x20ll) ? '1' : '0'), \
    (((i) & 0x10ll) ? '1' : '0'), \
    (((i) & 0x08ll) ? '1' : '0'), \
    (((i) & 0x04ll) ? '1' : '0'), \
    (((i) & 0x02ll) ? '1' : '0'), \
    (((i) & 0x01ll) ? '1' : '0')
#define PBI16(i) PBI08((i) >>  8), PBI08(i)
#define PBI32(i) PBI16((i) >> 16), PBI16(i)

#define PBI64(i) PBI32((i) >> 32), PBI32(i)
#define PBI(i) PBI64(i)

// end printf binary macros

#endif // WORD_H
