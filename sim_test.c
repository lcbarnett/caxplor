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
//	CLAP_CARG(n,       size_t,  10,           "length");
	CLAP_CARG(seed,    ulong,   0,            "random seed (or 0 for unpredictable)");
	puts("---------------------------------------------------------------------------------------\n");

	const size_t n = WBITS;

	double* const costab = calloc(n*n,sizeof(double));
	double* const sintab = calloc(n*n,sizeof(double));

	ft_tab(n,costab,sintab);

/*
	for (size_t i=0; i<n; ++i) {
		for (size_t j=0; j<n; ++j) printf("    % 8.6f",sintab[n*i+j]);
		putchar('\n');
	}
	putchar('\n');
*/
	for (size_t i=0; i<n; ++i) {
		for (size_t j=0; j<n; ++j) {
			if (costab[n*i+j] != costab[n*j+i]) printf("eek! i = %zu, j = %zu\n",i,j);
		}
	}

	mt_t rng;
	mt_seed(&rng,seed);
	const word_t w = wd_randomise(&rng);

	double* const wdftre = calloc(n,sizeof(double));
	double* const wdftim = calloc(n,sizeof(double));

	wd_dft(w,wdftre,wdftim,costab,sintab);

	const double nfac = 1.0/(double)n;
	for (size_t i=0; i<n; ++i) printf("    % 8.6f    % 8.6f\n",nfac*wdftre[i],nfac*wdftim[i]);
	putchar('\n');

	double* const wdspec = calloc(n,sizeof(double));
	for (size_t i=0; i<n; ++i) wdspec[i] = wdftre[i]*wdftre[i] + wdftim[i]*wdftim[i];

	for (size_t i=0; i<n; ++i) printf("    %8.6f\n",nfac*nfac*wdspec[i]);
	putchar('\n');

	const double ffac = (2.0*M_PI)/(double)n;
	FILE* gp = gp_popen(NULL,NULL);
	fprintf(gp,"set title \"Spatial Discrete Fourier Transform\"\n");
	fprintf(gp,"set xlabel \"frequency\"\n");
	fprintf(gp,"set ylabel \"DFT\"\n");
	fprintf(gp,"set xr [0.0:%g]\n",2.0*M_PI);
	fprintf(gp,"set yr [-1:1]\n");
	fprintf(gp,"set grid lt -1 lc \"grey\"\n");
	fprintf(gp,"plot \"-\" using 1:2 w lines t \"real\",  \"-\" using 1:2 w lines t \"imag\"\n");
	for (size_t i=0; i<n; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*wdftre[i]);
	fprintf(gp,"e\n");
	for (size_t i=0; i<n; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*wdftim[i]);
	fprintf(gp,"e\n");
	gp_pclose(gp);

	free(wdspec);
	free(wdftim);
	free(wdftre);
	free(sintab);
	free(costab);

	return EXIT_SUCCESS;
}
