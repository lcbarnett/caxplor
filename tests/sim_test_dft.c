#include "word.h"
#include "clap.h"
#include "utils.h"

void mw_dft_ref(const size_t n, const word_t* const w, double* const dftre, double* const dftim, double* const dps, const double* const costab)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	const double* const sintab = costab+m*m; // sin table is offset by m*m from cos table!
	for (size_t k=0;k<q;++k) {
		const size_t mk = m*k;
		double dftrek = 0.0;
		double dftimk = 0.0;
		for (size_t j=0;j<m;++j) {
			if (BITON(w[j/WBITS],j%WBITS)) {
				dftrek += costab[mk+j];
				dftimk += sintab[mk+j];
			}
		}
		dftre[k] = +dftrek;
		dftim[k] = -dftimk;
	}
	if (dps != NULL) sqmag(q,dps,dftre,dftim);  // calculate discrete power spectrum
}

void mw_autocov_ref(const size_t n, const word_t* const w, double* const ac)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t k=0;k<q;++k) {
		word_t ack = 0;
		for (size_t j=0;j<m;++j) {
			const size_t i = j+k < m ? j+k : j+k-m; // wrap!
			ack += BITON(w[i/WBITS],i%WBITS)&BITON(w[j/WBITS],j%WBITS);
		}
		ac[k] = (double)ack;
	}
}

int sim_test(int argc, char* argv[], int info)
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

	if (info) return EXIT_SUCCESS; // display switches and return

	mt_t rng;
	mt_seed(&rng,seed);

	const size_t m = WBITS*n;
	const size_t q = m/2+1; // fine, because WBITS even!

	double* const costab = dft_cstab_alloc(m); // note: sin table is offset by m*m from cos table!

	word_t* const w = mw_alloc(n);
	mw_randomiseb(n,w,wbias,&rng);
	mw_print(n,w);
	putchar('\n');
	putchar('\n');

	double* const dftre1 = calloc(q,sizeof(double));
	double* const dftim1 = calloc(q,sizeof(double));
	double* const ftdps1 = calloc(q,sizeof(double));

	double* const dftre  = calloc(q,sizeof(double));
	double* const dftim  = calloc(q,sizeof(double));
	double* const ftdps  = calloc(q,sizeof(double));

	double* const acov   = calloc(q,sizeof(double));
	double* const acdps  = calloc(q,sizeof(double));

	double* const acov1  = calloc(q,sizeof(double));
	double* const acdps1 = calloc(q,sizeof(double));

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

	printf("\nft-ft max abs diff = %.4e\n", maxabdiff(q,ftdps1,ftdps));
	printf("ac-ac max abs diff = %.4e\n\n", maxabdiff(q,acov1, acov ));
	printf("ft-ac max abs diff = %.4e\n\n", maxabdiff(q,ftdps, acdps));

//for (size_t k=0;k<q;++k) printf("%4zu  % 16.4f  % 16.4f\n",k,ftdps[k],acdps[k]);

	if (logs) {
		for (size_t i=0;i<q;++i) ftdps1[i] = ftdps1[i] > 0.0 ? log(ftdps1[i]) : NAN;
		for (size_t i=0;i<q;++i) ftdps [i] = ftdps [i] > 0.0 ? log(ftdps [i]) : NAN;
		for (size_t i=0;i<q;++i) acdps [i] = acdps [i] > 0.0 ? log(acdps [i]) : NAN;
	}

	if (disp) {
		const double nfac = logs ? 1.0 : 1.0/(double)m;
		const double ffac = M_PI/(double)(q-1);
		FILE* gp = gp_popen(NULL,"CA Xplorer Test",150,0);
		fprintf(gp,"set xr [0.0:%g]\n",M_PI);
		fprintf(gp,"set xtics ('0' 0.0,'{/Symbol p}/2' %g,'{/Symbol p}' %g)\n",M_PI/2.0,M_PI);
		fprintf(gp,"set grid lt -1 lc \"grey\"\n");
		fprintf(gp,"set multiplot layout 2,1\n");
		fprintf(gp,"set title \"Spatial Discrete Power Spectrum\"\n");
		fprintf(gp,"set xlabel \"frequency\"\n");
		fprintf(gp,"set ylabel \"DPS\"\n");
		if (logs) fprintf(gp,"set logs y\n"); else fprintf(gp,"set yr [0:*]\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines t \"ft\",  \"-\" using 1:2 w lines t \"ac\"\n");
		for (size_t i=1; i<q; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*ftdps[i]);
//		fprintf(gp,"%g %g\n",ffac*(double)m,nfac*ftdps[0]);
		fprintf(gp,"e\n");
		for (size_t i=1; i<q; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*acdps[i]);
//		fprintf(gp,"%g %g\n",ffac*(double)m,nfac*acdps[0]);
		fprintf(gp,"e\n");
		fprintf(gp,"set title \"Spatial autocorrelation\"\n");
		fprintf(gp,"set xlabel \"correlation length\"\n");
		fprintf(gp,"set ylabel \"autocorrelation\"\n");
		fprintf(gp,"set yr [*:*]\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines not\n");
		for (size_t i=1; i<q; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,acov[i]/acov[0]);
//		fprintf(gp,"%g %g\n",ffac*(double)m,1.0);
		fprintf(gp,"e\n");
		fprintf(gp,"unset multiplot\n");
		gp_pclose(gp);
	}

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
