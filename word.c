#include <math.h>

#include "word.h"
#include "utils.h"

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

void wd_dft(const word_t w, dft_float_t* const wdftre, dft_float_t* const wdftim, const dft_float_t* const costab)
{
	const dft_float_t* const sintab = costab+WBITS*WBITS; // sin table is offset by m*m from cos table!
	word_t wtest;
	for (size_t i=0;i<WBITS;++i) wdftre[i] = 0.0;
	wtest = w;
	for (int j=0;j<WBITS;++j) {
		if (WONE&(wtest>>=1)) { // test j-th bit
			const dft_float_t* const ctabj = costab+WBITS*j;
			for (int i=0;i<WBITS;++i) wdftre[i] += ctabj[i];
		}
	}
	for (size_t i=0;i<WBITS;++i) wdftim[i] = 0.0;
	wtest = w;
	for (int j=0;j<WBITS;++j) {
		if (WONE&(wtest>>=1)) { // test j-th bit
			const dft_float_t* const stabj = sintab+WBITS*j;
			for (int i=0;i<WBITS;++i) wdftim[i] += stabj[i];
		}
	}
}

void wd_autocov(const word_t w, dft_float_t* const wac)
{
	// note: reflect around WBITS/2
	for (int i=0;i<WBITS;++i) {
		word_t iwaci = 0;
		int j = 0;
		for (;j<WBITS-i;++j) iwaci += BITON(w,j+i      )&BITON(w,j);
		for (;j<WBITS  ;++j) iwaci += BITON(w,j+i-WBITS)&BITON(w,j); // wrap
		wac[i] = (dft_float_t)iwaci;
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

void mw_dft(const size_t n, const word_t* const w, dft_float_t* const dftre, dft_float_t* const dftim, dft_float_t* const dps, const dft_float_t* const costab)
{
	const size_t m = n*WBITS;
//	const size_t mh = m/2; // no problem: WBITS will always be even :-)
	const dft_float_t* const sintab = costab+m*m; // sin table is offset by m*m from cos table!
	for (size_t i=0;i<m;++i) dftre[i] = (dft_float_t)0;
	for (size_t i=0;i<m;++i) dftim[i] = (dft_float_t)0;
	for (size_t jj=0,jdx=0;jj<n;++jj) {
		word_t wjj = w[jj];
		for (int j=0;j<WBITS;++j,++jdx,wjj>>=1) {
			if (WONE&wjj) {
				for (size_t i=0,ij=m*jdx;i<m;++i,++ij) {
					dftre[i] += costab[ij];
					dftim[i] -= sintab[ij];
				}
			}
		}
	}
#ifdef DFT_SINGLE_PREC_FLOAT
	if (dps != NULL) sqmagf(m,dps,dftre,dftim); // discrete power spectrum
#else
	if (dps != NULL) sqmag(m,dps,dftre,dftim);  // discrete power spectrum
#endif
}

void mw_autocov(const size_t n, const word_t* const w, dft_float_t* const ac)
{
	size_t idx = 0;
	for (size_t ii=0;ii<n;++ii) {
		for (int i=0;i<WBITS;++i,++idx) {
			word_t aci = 0;
			size_t jii = ii;
			int    ji  = i;
			for (size_t jj=0;jj<n;++jj) {
				for (int j=0;j<WBITS;++j) {
					aci += BITON(w[jj],j)&BITON(w[jii],ji);
					if (++ji == WBITS) {
						if (++jii == n) jii = 0; // wrap word
						ji = 0;
					}
				}
			}
			ac[idx] = (dft_float_t)aci;
		}
	}
}

void mw_automi(const size_t n, const word_t* const w, double* const ami)
{
	const size_t m = n*WBITS;
	const double fac = 1.0/(double)m;
	const double p0 = fac*mw_nsetbits(n,w);
	const double entro = -xlog2x(p0)-xlog2x(1.0-p0);
	int bin[4];
	ami[0] = entro;
	for (size_t k=1;k<m;++k) {
		for (size_t r=0;r<4;++r) bin[r] = 0;
		for (size_t j=0;j<m;++j) {
			const size_t i = j+k < m ? j+k : j+k-m; // wrap!
			++bin[MIIDX(w[i/WBITS],i%WBITS,w[j/WBITS],j%WBITS)];
		}
		ami[k] = 2.0*entro+xlog2x(fac*(double)bin[0])+xlog2x(fac*(double)bin[1])+xlog2x(fac*(double)bin[2])+xlog2x(fac*(double)bin[3]);
	}
}

/*********************************************************************/
