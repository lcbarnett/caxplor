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
	CLAP_CARG(wbias,   double,  0.5,          "random word bias");
	CLAP_CARG(seed,    ulong,   0,            "random seed (or 0 for unpredictable)");
	CLAP_CARG(disp,    int,     1,            "display result?");
	puts("---------------------------------------------------------------------------------------\n");

	mt_t rng;
	mt_seed(&rng,seed);

	const size_t m = WBITS*n;

	double* const costab = dft_costab_alloc(m);
	double* const sintab = dft_sintab_alloc(m);

	word_t* const w = mw_alloc(n);
	mw_randomiseb(n,w,wbias,&rng);
	mw_print(n,w);
	putchar('\n');
	putchar('\n');

	double* const dftre = calloc(m,sizeof(double));
	double* const dftim = calloc(m,sizeof(double));
	double* const ftdps = calloc(m,sizeof(double));

	double* const acov  = calloc(m,sizeof(double));
	double* const acdps = calloc(m,sizeof(double));

	double ts,te;

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	mw_dft(n,w,dftre,dftim,costab,sintab);
	for (size_t i=0; i<m; ++i) ftdps[i] = dftre[i]*dftre[i] + dftim[i]*dftim[i];
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("ft dps time = %8.6f\n",te-ts);

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	mw_autocov(n,w,acov);
	dps(m,acdps,acov,costab);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("ac dps time = %8.6f\n",te-ts);

//	for (size_t i=0; i<m; ++i) printf("%4zu   %4.0f    % 12.4f\n",i,acov[i],acdps[i]);

	if (disp) {
		const double nfac = 1.0/(double)m;
		const double ffac = (2.0*M_PI)/(double)m;
		FILE* gp = gp_popen(NULL,NULL);
		fprintf(gp,"set title \"Spatial Discrete Power Spectrum\"\n");
		fprintf(gp,"set xlabel \"frequency\"\n");
		fprintf(gp,"set ylabel \"DPS\"\n");
		fprintf(gp,"set xr [0.0:%g]\n",2.0*M_PI);
		fprintf(gp,"# set yr [-1:1]\n");
		fprintf(gp,"set grid lt -1 lc \"grey\"\n");
		fprintf(gp,"plot \"-\" using 1:2 w lines t \"ft\",  \"-\" using 1:2 w lines t \"ac\"\n");
		for (size_t i=0; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*ftdps[i]);
		fprintf(gp,"e\n");
		for (size_t i=0; i<m; ++i) fprintf(gp,"%g %g\n",ffac*(double)i,nfac*acdps[i]);
		fprintf(gp,"e\n");
		gp_pclose(gp);
	}

	free(acdps);
	free(acov);
	free(ftdps);
	free(dftim);
	free(dftre);
	free(w);
	free(sintab);
	free(costab);
	return EXIT_SUCCESS;
}
