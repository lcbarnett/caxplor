#include <math.h>

#include "word.h"
#include "utils.h"

/*********************************************************************/
/*                      miscellaneous                                */
/*********************************************************************/

void ft_tab(const size_t n, double* const costab, double* const sintab)
{
	const double fac = (2.0*M_PI)/(double)n;
	for (size_t i=0; i<n; ++i) {
		double* const ctni = costab+n*i;
		double* const stni = sintab+n*i;
		for (size_t j=i; j<n; ++j) sincos(fac*(double)(i*j),stni+j,ctni+j);
	}
	// symmetrise across diagonal
	for (size_t i=0; i<n; ++i) {
		double* const cti  = costab+i;
		double* const ctni = costab+n*i;
		for (size_t j=i+1; j<n; ++j) *(cti+n*j) = *(ctni+j);
	}
	for (size_t i=0; i<n; ++i) {
		double* const sti  = sintab+i;
		double* const stni = sintab+n*i;
		for (size_t j=i+1; j<n; ++j) *(sti+n*j) = *(stni+j);
	}
}

/*********************************************************************/
/*                      single-word                                  */
/*********************************************************************/

void wd_print(const word_t w)
{
	for (word_t m=WHIONE;m!=0;m>>=1) putchar(NZWORD(w&m)?'1':'0');
}

void wd_fprint(const word_t w, FILE* const fstream)
{
	for (word_t m=WHIONE;m!=0;m>>=1) fputc(NZWORD(w&m)?'1':'0',fstream);
}

void wd_prints(const word_t w)
{
	for (word_t m=WHIONE;m!=0;m>>=1) {putchar(' '); putchar(NZWORD(w&m)?'1':'0');}
}

void wd_fprints(const word_t w, FILE* const fstream)
{
	for (word_t m=WHIONE;m!=0;m>>=1) {fputc(' ',fstream); fputc(NZWORD(w&m)?'1':'0',fstream);}
}

void wd_printc(const word_t w, const int W)
{
	for (int i=0;i<WBITS;++i) {
		putchar(NZWORD((w<<i)&WHIONE)?'1':'0');
		if (!((i+1)%W) && i < WBITS1) putchar(',');
	}
}

void wd_fprintc(const word_t w, const int W, FILE* const fstream)
{
	for (int i=0;i<WBITS;++i) {
		fputc(NZWORD((w<<i)&WHIONE)?'1':'0'? '1':'0',fstream);
		if (!((i+1)%W) && i < WBITS1) fputc(',',fstream);
	}
}

void wd_print_lo(const word_t w, const int b)
{
	for (word_t m=(WONE<<(b-1));m!=0;m>>=1) putchar(NZWORD(w&m)?'1':'0');
}

void wd_fprint_lo(const word_t w, const int b, FILE* const fstream)
{
	for (word_t m=(WONE<<(b-1));m!=0;m>>=1) fputc(NZWORD(w&m)?'1':'0',fstream);
}

void wd_dft(const word_t w, double* const wdftre, double* const wdftim, const double* const costab, const double* const sintab)
{
	word_t wtest;
	for (size_t i=0;i<WBITS;++i) wdftre[i] = 0.0;
	wtest = w;
	for (int j=0;j<WBITS;++j) {
		if (WONE&(wtest>>=1)) { // test j-th bit
			const double* const ctabj = costab+WBITS*j;
			for (int i=0;i<WBITS;++i) wdftre[i] += ctabj[i];
		}
	}
	for (size_t i=0;i<WBITS;++i) wdftim[i] = 0.0;
	wtest = w;
	for (int j=0;j<WBITS;++j) {
		if (WONE&(wtest>>=1)) { // test j-th bit
			const double* const stabj = sintab+WBITS*j;
			for (int i=0;i<WBITS;++i) wdftim[i] += stabj[i];
		}
	}
}

/*********************************************************************/
/*                      multi-word                                   */
/*********************************************************************/

word_t* mw_alloc(const size_t n)
{
	word_t* const w =  calloc(n,sizeof(word_t)); // note: calloc zero-initialises - this is a Good Thing
	PASSERT(w != NULL,"memory allocation failed");
	return w;

}

word_t* mw_alloc_copy(const size_t n, const word_t* const wsrc)
{
	word_t* const wdest = malloc(n*sizeof(word_t));
	PASSERT(wdest != NULL,"memory allocation failed");
	memcpy(wdest,wsrc,n*sizeof(word_t));
	return wdest;
}

void mw_fprint(const size_t n, const word_t* const w, FILE* const fstream)
{
	for (const word_t* u=w+n-1;u>=w;--u) wd_fprint(*u,fstream);
}

void mw_print(const size_t n, const word_t* const w)
{
	mw_fprint(n,w,stdout);
}

void mw_fprints(const size_t n, const word_t* const w, FILE* const fstream)
{
	for (const word_t* u=w+n-1;u>=w;--u) wd_fprints(*u,fstream);
}

void mw_prints(const size_t n, const word_t* const w)
{
	mw_fprints(n,w,stdout);
}

void mw_fprint_bin(const size_t n, const word_t* const w, FILE* const fstream)
{
	const int bchunk = 32;
	const int WB = WBITS-bchunk;
	int firstword = 1;
	for (const word_t* u=w+n-1;u>=w;--u) {
		if (firstword) {
			printf("%"PRIw,(*u)>>WB);
			for (int i=WBITS-bchunk;i!=0;i-=bchunk) {
				printf(",%"PRIw,((*u)<<(WBITS-i))>>WB);
			}
			firstword = 0;
		}
		else {
			for (int i=WBITS;i!=0;i-=bchunk) {
				printf(",%"PRIw,((*u)<<(WBITS-i))>>WB);
			}
		}
	}
}

void mw_print_bin(const size_t n, const word_t* const w)
{
	mw_fprint_bin(n,w,stdout);
}

void mw_run(const size_t I, const size_t n, word_t* const w, const int B, const word_t* const f)
{
	if (I == 0) return; // do nothing
	const size_t J = I/2;
	word_t ww[n];
	for (size_t j=0;j<J;++j) {
		mw_filter(n,ww,w,B,f);
		mw_filter(n,w,ww,B,f);
	}
	if (2*J != I) { // odd number of iterations - one more to go
		mw_filter(n,ww,w,B,f);
		mw_copy(n,w,ww);
	}
}

void mw_dft(const size_t n, const word_t* const w, double* const wdftre, double* const wdftim, const double* const costab, const double* const sintab)
{
	const size_t m = n*WBITS;
	for (size_t i=0;i<m;++i) wdftre[i] = 0.0;
	for (size_t i=0;i<m;++i) wdftim[i] = 0.0;
	for (size_t l=0;l<n;++l) {
		word_t wl = w[l];
		for (int j=0;j<WBITS;++j) {
			if (WONE&(wl>>=1)) {
				const double* const costabj = costab+(l*WBITS+(size_t)j)*m;
				const double* const sintabj = sintab+(l*WBITS+(size_t)j)*m;
				for (size_t k=0;k<n;++k) {
					const double* const costabjk = costabj+k*WBITS;
					const double* const sintabjk = sintabj+k*WBITS;
					double* const wdftrek = wdftre+k*WBITS;
					double* const wdftimk = wdftim+k*WBITS;
					for (int i=0;i<WBITS;++i) wdftrek[i] += costabjk[i];
					for (int i=0;i<WBITS;++i) wdftimk[i] -= sintabjk[i];
				}
			}
		}
	}
}

void mw_dft_ref(const size_t n, const word_t* const w, double* const wdftre, double* const wdftim, const double* const costab, const double* const sintab)
{
	const size_t m = n*WBITS;
	for (size_t i=0;i<m;++i) wdftre[i] = 0.0;
	for (size_t i=0;i<m;++i) wdftim[i] = 0.0;
	for (size_t l=0;l<n;++l) {
		word_t wl = w[l];
		for (int j=0;j<WBITS;++j) {
			if (WONE&(wl>>=1)) {
				for (size_t k=0;k<n;++k) {
					for (int i=0;i<WBITS;++i) wdftre[k*WBITS+(size_t)i] += costab[(l*WBITS+(size_t)j)*m+k*WBITS+(size_t)i];
					for (int i=0;i<WBITS;++i) wdftim[k*WBITS+(size_t)i] -= sintab[(l*WBITS+(size_t)j)*m+k*WBITS+(size_t)i];
				}
			}
		}
	}
}

/*********************************************************************/
