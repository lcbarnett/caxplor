#include "word.h"
#include "clap.h"
#include "utils.h"

void mw_dft_ref(const size_t n, const word_t* const w, dft_float_t* const dftre, dft_float_t* const dftim, dft_float_t* const dps, const dft_float_t* const costab)
{
	const size_t m = n*WBITS;
	const dft_float_t* const sintab = costab+m*m; // sin table is offset by m*m from cos table!
	for (size_t idx=0,ii=0;ii<n;++ii) {
		for (int i=0;i<WBITS;++i,++idx) {
			const size_t midx = m*idx;
			dft_float_t dftrei = 0.0;
			dft_float_t dftimi = 0.0;
			for (size_t jdx=0,jj=0;jj<n;++jj) {
				for (int j=0;j<WBITS;++j,++jdx) {
					dftrei += BITON(w[jj],j) ? (dft_float_t)costab[midx+jdx] : (dft_float_t)0;
					dftimi += BITON(w[jj],j) ? (dft_float_t)sintab[midx+jdx] : (dft_float_t)0;
				}
			}
			dftre[idx] = +dftrei;
			dftim[idx] = -dftimi;
		}
	}
#ifdef DFT_SINGLE_PREC_FLOAT
	if (dps != NULL) sqmagf(m,dps,dftre,dftim); // calculate discrete power spectrum
#else
	if (dps != NULL) sqmag(m,dps,dftre,dftim);  // calculate discrete power spectrum
#endif
}

void mw_autocov_ref(const size_t n, const word_t* const w, dft_float_t* const ac)
{
	const size_t m = n*WBITS;
	for (size_t k=0;k<m;++k) {
		word_t ack = (dft_float_t)0;
		for (size_t j=0;j<m;++j) {
			const size_t i   = j+k < m ? j+k : j+k-m; // wrap!
			const size_t jj  = j/WBITS;
			const int    jjj = j%WBITS;
			const size_t ii  = i/WBITS;
			const int    iii = i%WBITS;
			ack += BITON(w[jj],jjj)&BITON(w[ii],iii);
		}
		ac[k] = (dft_float_t)ack;
	}
}

void mw_dps(const size_t n, const word_t* const w, dft_float_t* const dps, const dft_float_t* const costab)
{
	//NOTE: not very efficient - use mw_dft instead!
	const size_t m = n*WBITS;
	for (size_t k=0;k<m;++k) dps[k] = (dft_float_t)0;
	for (size_t i=0,ii=0;ii<n;++ii) {
		for (int iii=0;iii<WBITS;++iii,++i) {
			if (BITON(w[ii],i)) {
				for (size_t j=0,jj=0;jj<n;++jj) {
					for (int jjj=0;jjj<WBITS;++jjj,++j) {
						if (BITON(w[jj],j)) {
							const size_t md = m*(i>j?i-j:j-i);
							for (size_t k=0;k<m;++k) dps[k] += costab[md+k];
						}
					}
				}
			}
		}
	}
}

void mw_dps_ref(const size_t n, const word_t* const w, dft_float_t* const dps, const dft_float_t* const costab)
{
	const size_t m = n*WBITS;
	for (size_t k=0,kk=0;kk<n;++kk) {
		for (int kkk=0;kkk<WBITS;++kkk,++k) {
			const size_t mk = m*k;
			dft_float_t dpsk = (dft_float_t)0;
			for (size_t i=0,ii=0;ii<n;++ii) {
				for (int iii=0;iii<WBITS;++iii,++i) {
					if (BITON(w[ii],i)) {
						for (size_t j=0,jj=0;jj<n;++jj) {
							for (int jjj=0;jjj<WBITS;++jjj,++j) {
								const size_t d = i>j ? i-j : j-i;
								if (BITON(w[jj],j)) dpsk += costab[mk+d];
							}
						}
					}
				}
			}
			dps[k] = dpsk;
		}
	}
}

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(n,       size_t,  3,            "number of words");
	CLAP_CARG(wbias,   double,  0.5,          "random word bias");
	CLAP_CARG(seed,    ulong,   0,            "random seed (or 0 for unpredictable)");
	CLAP_CARG(disp,    int,     1,            "display result?");
	CLAP_CARG(logs,    int,     0,            "log scale?");
	puts("---------------------------------------------------------------------------------------\n");

	mt_t rng;
	mt_seed(&rng,seed);

	const size_t m = WBITS*n;

	dft_float_t* const costab = dft_cstab_alloc(m); // note: sin table is offset by m*m from cos table!

	word_t* const w = mw_alloc(n);
	mw_randomiseb(n,w,wbias,&rng);
	mw_print(n,w);
	putchar('\n');
	putchar('\n');

	dft_float_t* const dftre1 = calloc(m,sizeof(dft_float_t));
	dft_float_t* const dftim1 = calloc(m,sizeof(dft_float_t));
	dft_float_t* const ftdps1 = calloc(m,sizeof(dft_float_t));

	dft_float_t* const dftre  = calloc(m,sizeof(dft_float_t));
	dft_float_t* const dftim  = calloc(m,sizeof(dft_float_t));
	dft_float_t* const ftdps  = calloc(m,sizeof(dft_float_t));

	dft_float_t* const dcdps  = calloc(m,sizeof(dft_float_t));

	dft_float_t* const acov   = calloc(m,sizeof(dft_float_t));
	dft_float_t* const acdps  = calloc(m,sizeof(dft_float_t));

	dft_float_t* const acov1  = calloc(m,sizeof(dft_float_t));
	dft_float_t* const acdps1 = calloc(m,sizeof(dft_float_t));

	double ts,te;

	ts = timer();
	mw_dft_ref(n,w,dftre1,dftim1,ftdps1,costab);
	te = timer();
	printf("ft ref time = %8.6f\n",te-ts);

	ts = timer();
	mw_dft(n,w,dftre,dftim,ftdps,costab);
	te = timer();
	printf("ft dps time = %8.6f\n",te-ts);

	ts = timer();
	mw_autocov_ref(n,w,acov1);
	ac2dps(m,acdps1,acov1,costab);
	te = timer();
	printf("\nac ref time = %8.6f\n",te-ts);

	ts = timer();
	mw_autocov(n,w,acov);
	ac2dps(m,acdps,acov,costab);
	te = timer();
	printf("ac dps time = %8.6f\n\n",te-ts);

	ts = timer();
	mw_dps(n,w,dcdps,costab);
	te = timer();
	printf("dc dps time = %8.6f\n",te-ts);


	printf("\nft-ft max abs diff = %.4e\n", maxabdifff(m,ftdps1,ftdps));
	printf("ac-ac max abs diff = %.4e\n\n", maxabdifff(m,acov1,acov));

	printf("ft-ac max abs diff = %.4e\n",   maxabdifff(m,ftdps, acdps));
	printf("ft-dc max abs diff = %.4e\n\n", maxabdifff(m,ftdps, dcdps));


//const int nsb = mw_nsetbits(n,w);
//printf("\nset bits = %d\n\n",nsb);
//for (size_t i=0;i<m;++i) printf("%4zu  % 12.4f  % 12.4f\n",i,acov[i],acov1[i]);


// for (size_t i=0;i<m;++i) printf("%4zu  % 12.4f  % 12.4f\n",i,ftdps[i],dcdps[i]);
// for (size_t i=0; i<m; ++i) printf("%4zu  % 12.4f  % 12.4f\n",i,dftre[i]-dftre1[i],dftim[i]-dftim1[i]);


	if (logs) {
		for (size_t i=0;i<m;++i) ftdps1[i] = ftdps1[i] > 0.0 ? logf(ftdps1[i]) : NAN;
		for (size_t i=0;i<m;++i) ftdps [i] = ftdps [i] > 0.0 ? logf(ftdps [i]) : NAN;
		for (size_t i=0;i<m;++i) acdps [i] = acdps [i] > 0.0 ? logf(acdps [i]) : NAN;
		for (size_t i=0;i<m;++i) dcdps [i] = dcdps [i] > 0.0 ? logf(dcdps [i]) : NAN;
	}



	if (disp) {
		const dft_float_t nfac = logs ? (dft_float_t)1 : (dft_float_t)1/(dft_float_t)m;
		const dft_float_t ffac = (dft_float_t)(2.0*M_PI)/(dft_float_t)m;
		FILE* gp = gp_popen(NULL,"wxt size 640,960 nobackground enhanced title 'CA Explorer Test' persist raise");
		fprintf(gp,"set xr [0.0:%g]\n",2.0*M_PI);
		fprintf(gp,"set xtics ('0' 0.0,'{/Symbol p}/2' %g,'{/Symbol p}' %g,'3{/Symbol p}/2' %g,'2{/Symbol p}' %g)\n",M_PI/2.0,M_PI,3.0*M_PI/2.0,2*M_PI);
		fprintf(gp,"set grid lt -1 lc \"grey\"\n");
		fprintf(gp,"set multiplot layout 2,1\n");
		fprintf(gp,"set title \"Spatial Discrete Power Spectrum\"\n");
		fprintf(gp,"set xlabel \"frequency\"\n");
		fprintf(gp,"set ylabel \"DPS\"\n");
		if (logs) fprintf(gp,"set logs y\n"); else fprintf(gp,"set yr [0:*]\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines t \"ft\",  \"-\" using 1:2 w lines t \"ac\"\n");
		for (size_t i=1; i<m-1; ++i) fprintf(gp,"%g %g\n",ffac*(dft_float_t)i,nfac*ftdps[i]);
//		fprintf(gp,"%g %g\n",ffac*(double)m,nfac*ftdps[0]);
		fprintf(gp,"e\n");
		for (size_t i=1; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(dft_float_t)i,nfac*acdps[i]);
//		fprintf(gp,"%g %g\n",ffac*(double)m,nfac*acdps[0]);
		fprintf(gp,"e\n");
		fprintf(gp,"set title \"Spatial autocorrelation\"\n");
		fprintf(gp,"set xlabel \"correlation length\"\n");
		fprintf(gp,"set ylabel \"autocorrelation\"\n");
		fprintf(gp,"set yr [*:*]\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines not\n");
		for (size_t i=1; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(dft_float_t)i,acov[i]/acov[0]);
//		fprintf(gp,"%g %g\n",ffac*(dft_float_t)m,1.0);
		fprintf(gp,"e\n");
		fprintf(gp,"unset multiplot\n");
		gp_pclose(gp);
	}

	free(dcdps);
	free(acdps);
	free(acdps1);
	free(acov);
	free(acov1);
	free(ftdps);
	free(dftim);
	free(dftre);
	free(ftdps1);
	free(dftim1);
	free(dftre1);
	free(w);
	free(costab);
	return EXIT_SUCCESS;
}
