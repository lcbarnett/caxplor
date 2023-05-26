#include <time.h>

#include "ca.h"
#include "rtab.h"
#include "clap.h"
#include "strman.h"

int sim_ana(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(rsiz,    int,     5,            "CA rule size");
	CLAP_VARG(bmax,    size_t,  0,            "CA max bits to set (or zero for half)");
	CLAP_CARG(bseed,   ulong,   0,            "CA rule random seed (0 for unpredictable)");
	CLAP_CARG(m,       int,     20,           "CA sequence length for entropy calculation");
	CLAP_CARG(iters,   int,     1,            "CA iterations before entropy calculation");
	CLAP_CARG(S,       size_t,  100,          "sample size");
	puts("---------------------------------------------------------------------------------------\n");

	sm_create(sm);

	// pseudo-random number generators
	mt_t brng;
	mt_seed(&brng,bseed);

	// allocate CA rule table
	word_t* const rtab = rt_alloc(rsiz);

	// allocate storage
	bmax = (bmax == 0 ? POW2(rsiz-1) : bmax);
	double*  const HH = calloc((bmax+1)*S,sizeof(double));
	double** const H  = calloc(bmax+1,sizeof(double*));
	for (size_t b=0; b<=bmax; ++b) H[b] = HH+b*S;

	// entropy calculation
	const double oom = 1.0/(double)m;
	for (size_t b=0; b<=bmax; ++b) {
		const double ts = (double)clock()/(double)CLOCKS_PER_SEC;
		printf("b = %3zu of %zu ...",b,bmax); fflush(stdout);
		for (size_t s=0; s<S; ++s) {
			rt_randomb(rsiz,rtab,b,&brng);
			H[b][s] = oom*rt_entro(rsiz,rtab,m,iters);
		}
		const double et = (double)clock()/(double)CLOCKS_PER_SEC - ts;
		const int mins = (int)floor(et/60.0);
		const double secs = fmod(et,60.0);
		printf(" %3d mins, %5.2f secs\n",mins,secs);
	}

	char* const gpfile = sm_printf(sm,"caana_B%d_m%d_i%d_S%zu",rsiz,m,iters,S);
	FILE* const gp = gp_fopen(gpfile,NULL,NULL);
	fprintf(gp,"set title 'CA entropy vs bias (rsiz = %d, m = %d, iters = %d, samples = %zu)'\n",rsiz,m,iters,S);
	fprintf(gp,"set xlabel 'rule bits set'\n");
	fprintf(gp,"set ylabel 'normalised entropy'\n");
	fprintf(gp,"set xr[-0.5:%g]\n",(double)bmax+0.5);
	fprintf(gp,"set yr[0.0:1.0]\n");
	fprintf(gp,"set ytics 0.1\n");
	fprintf(gp,"set grid\n");
	fprintf(gp,"plot '-' u 1:2 w points pt 7 ps 0.75 not\n");
	for (size_t b=0; b<=bmax; ++b) {
		for (size_t s=0; s<S; ++s) {
			fprintf(gp,"%zu\t%g\n",b,H[b][s]);
		}
	}
	fprintf(gp,"e\n");
	if (fclose(gp) == -1) PEEXIT("failed to close Gnuplot data file\n");
	gp_fplot(gpfile,NULL,NULL);

	// free storage
	free(H);
	free(HH);
	free(rtab);

	sm_destroy(sm);

	return EXIT_SUCCESS;
}
