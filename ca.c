#include "ca.h"
#include "utils.h"

/*********************************************************************/
/*                      CA (multi-word)                              */
/*********************************************************************/

void ca_fprint(const size_t I, const size_t n, const word_t* const ca, FILE* const fstream)
{
	for (const word_t* w=ca;w<ca+I*n;w+=n) {mw_fprint(n,w,fstream); fnewline(fstream);}
}

void ca_print(const size_t I, const size_t n, const word_t* const ca)
{
	ca_fprint(I,n,ca,stdout);
}

void ca_fprints(const size_t I, const size_t n, const word_t* const ca, FILE* const fstream)
{
	for (const word_t* w=ca;w<ca+I*n;w+=n) {mw_fprints(n,w,fstream); fnewline(fstream);}
}

void ca_prints(const size_t I, const size_t n, const word_t* const ca)
{
	ca_fprints(I,n,ca,stdout);
}

void ca_run(const size_t I, const size_t n, word_t* const ca, word_t* const cawrk, const int B, const word_t* const rtab, const int uto)
{
	for (word_t* w=ca+n;w<ca+I*n;w+=n) mw_filter(n,w,w-n,B,rtab);
	if (uto) {
		ASSERT(cawrk != NULL,"Need CA work buffer to unwrap!");
		mw_copy(I*n,cawrk,ca);
		ca_rotl(I,n,ca,cawrk,uto);
	}
}

void ca_filter(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int B, const word_t* const rtab)
{
	word_t* wnew = ca;
	for (const word_t* w=caold;w<caold+I*n;w+=n,wnew+=n) mw_filter(n,wnew,w,B,rtab);
}

size_t ca_period(const size_t I, const size_t n, const word_t* const ca, const int B, const word_t* const rtab, int* const rot)
{
	word_t mword1[n];
	word_t mword2[n];
	word_t* wold = mword1;
	word_t* wnew = mword2;
	mw_copy(n,wold,ca);
	for (size_t i=0;i<I;++i) {
		mw_filter(n,wnew,wold,B,rtab);
		*rot = mw_equiv(n,ca,wnew);
		if (*rot >= 0) return i+1;
		SWAP(word_t*,wnew,wold);
	}
	return I;
}

void ca_rotl(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits)
{
	size_t nb = (size_t)nbits;
	size_t b = 0;
	for (size_t r=0;r<I*n;r+=n,b+=nb) mw_rotl(n,ca+r,caold+r,b);
}

void ca_rotr(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits)
{
	size_t nb = (size_t)nbits;
	size_t b = 0;
	for (size_t r=0; r<I*n; r+=n, b+=nb) mw_rotr(n,ca+r,caold+r,b);
}

void ca_reverse(const size_t I, const size_t n, word_t* const ca, const word_t* const caold)
{
	for (size_t r=0; r<I*n; r += n) mw_reverse(n,ca+r,caold+r);
}

void ca_part_count(const size_t I, const size_t n, const word_t* const ca, const int P, const size_t nparts, partf_t* const ppw, const int sortem)
{
	ASSERT(nparts == POW2(P),"number of particles must equal 2^(particle breadth)");
	word_t* const wcgrain = mw_alloc(n);
	const word_t* const caend = ca+I*n;
	for (word_t part=0;part<nparts;++part) {
		ulong ppwp = 0;
		for (const word_t* w=ca;w<caend;w+=n) {
			mw_zero(n,wcgrain);
			mw_parts(n,wcgrain,w,P,part);
			ppwp += (ulong)mw_nsetbits(n,wcgrain);
		}
		ppw[part].w = part;
		ppw[part].f = (double)ppwp/(double)I;
	}
	if (sortem) qsort(ppw,nparts,sizeof(partf_t),partf_comp);
	free(wcgrain);
}

void ca_dps(const size_t I, const size_t n, const word_t* const ca, double* const dps, const double* const costab)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	double* const dftre = calloc(q,sizeof(double));
	TEST_ALLOC(dftre);
	double* const dftim = calloc(q,sizeof(double));
	TEST_ALLOC(dftim);
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const dpsr = dps+q*row;
		mw_dft(n,car,dftre,dftim,dpsr,costab);
	}
	free(dftim);
	free(dftre);
}

void ca_autocov(const size_t I, const size_t n, const word_t* const ca, double* const ac)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const acr = ac+q*row;
		mw_autocov(n,car,acr);
	}
}

void ca_automi(const size_t I, const size_t n, const word_t* const ca, double* const ami)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const amir = ami+q*row;
		mw_automi(n,car,amir);
	}
}
