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
	puts("---------------------------------------------------------------------------------------\n");

	const size_t m = WBITS*n;

	double* const costab = calloc(m*m,sizeof(double));
	double* const sintab = calloc(m*m,sizeof(double));

	ft_tab(m,costab,sintab);

	for (size_t i=0; i<m; ++i) {
		for (size_t j=0; j<m; ++j) {
			if (costab[m*i+j] != costab[m*j+i]) printf("eek! i = %zu, j = %zu\n",i,j);
		}
	}

	mt_t rng;
	mt_seed(&rng,seed);
	word_t* const w = mw_alloc(n);
	mw_randomise(n,w,&rng);

	double* const wdftre = calloc(m,sizeof(double));
	double* const wdftim = calloc(m,sizeof(double));

	double* const wdftre1 = calloc(m,sizeof(double));
	double* const wdftim1 = calloc(m,sizeof(double));

	mw_dft    (n,w,wdftre, wdftim, costab,sintab);
	mw_dft_ref(n,w,wdftre1,wdftim1,costab,sintab);

	for (size_t i=0; i<m; ++i) if (wdftre1[i] != wdftre[i]) printf("real : i = %zu, diff = %e\n",i,fabs(wdftre1[i]-wdftre[i]));
	for (size_t i=0; i<m; ++i) if (wdftim1[i] != wdftim[i]) printf("imag : i = %zu, diff = %e\n",i,fabs(wdftim1[i]-wdftim[i]));


	double* const wdspec = calloc(m,sizeof(double));
	for (size_t i=0; i<m; ++i) wdspec[i] = wdftre[i]*wdftre[i] + wdftim[i]*wdftim[i];

	// for (size_t i=0; i<m; ++i) printf("    %8.6f\n",nfac*nfac*wdspec[i]);
	// putchar('\n');

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

	free(wdspec);
	free(wdftim1);
	free(wdftre1);
	free(wdftim);
	free(wdftre);
	free(w);
	free(sintab);
	free(costab);

	return EXIT_SUCCESS;
}
