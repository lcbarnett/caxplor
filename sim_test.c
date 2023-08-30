#include <time.h>

#include "word.h"
#include "clap.h"
#include "utils.h"

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(n,       size_t,  3,            "number of words");
	CLAP_CARG(seed,    ulong,   0,            "random seed (or 0 for unpredictable)");
	CLAP_CARG(disp,    int,     1,            "display result?");
	puts("---------------------------------------------------------------------------------------\n");

	mt_t rng;
	mt_seed(&rng,seed);
/*
	double* const wac1 = calloc(WBITS,sizeof(double));
	double* const wac2 = calloc(WBITS,sizeof(double));

	word_t w = wd_randomise(&rng);

	wd_autocov_ref(w,wac1);
	wd_autocov    (w,wac2);

	for (size_t i=0; i<WBITS; ++i) {
		printf("%12.0f  %12.0f",wac1[i],wac2[i]);
		if (wac1[i] != wac2[i]) printf("    ***");
		putchar('\n');
	}

	free(wac2);
	free(wac1);
*/
/*
	const size_t m = WBITS*n;

	double* const costab = dft_costab_alloc(m);
	double* const sintab = dft_sintab_alloc(m);

	word_t* const w = mw_alloc(n);
	mw_randomise(n,w,&rng);

	double* const wdftre = calloc(m,sizeof(double));
	double* const wdftim = calloc(m,sizeof(double));

	double* const wdftre1 = calloc(m,sizeof(double));
	double* const wdftim1 = calloc(m,sizeof(double));

	double ts,te;
	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	mw_dft_ref(n,w,wdftre1,wdftim1,costab,sintab);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("dft ref time = %8.6f\n",te-ts);
	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	mw_dft(n,w,wdftre, wdftim, costab,sintab);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("dft new time = %8.6f\n\n",te-ts);

	for (size_t i=0; i<m; ++i) if (wdftre1[i] != wdftre[i]) printf("real : i = %zu, diff = %e\n",i,fabs(wdftre1[i]-wdftre[i]));
	for (size_t i=0; i<m; ++i) if (wdftim1[i] != wdftim[i]) printf("imag : i = %zu, diff = %e\n",i,fabs(wdftim1[i]-wdftim[i]));

	double* const wdspec = calloc(m,sizeof(double));
	for (size_t i=0; i<m; ++i) wdspec[i] = wdftre[i]*wdftre[i] + wdftim[i]*wdftim[i];

	if (disp) {
		const double nfac = 1.0/(double)m;
		const double ffac = (2.0*M_PI)/(double)m;
		FILE* gp = gp_popen(NULL,NULL);
		fprintf(gp,"set title \"Spatial Discrete Fourier Transform\"\n");
		fprintf(gp,"set xlabel \"frequency\"\n");
		fprintf(gp,"set ylabel \"DFT\"\n");
		fprintf(gp,"set xr [0.0:%g]\n",2.0*M_PI);
		fprintf(gp,"set yr [-1:1]\n");
		fprintf(gp,"set grid lt -1 lc \"grey\"\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines t \"real\",  \"-\" using 1:2 w lines t \"imag\"\n");
		for (size_t i=0; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*wdftre[i]);
		fprintf(gp,"e\n");
		for (size_t i=0; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*wdftim[i]);
		fprintf(gp,"e\n");
		gp_pclose(gp);
	}

	free(wdspec);
	free(wdftim1);
	free(wdftre1);
	free(wdftim);
	free(wdftre);
	free(w);
	free(sintab);
	free(costab);
*/
	return EXIT_SUCCESS;
}
