#include <time.h>

#include "ca.h"
#include "rtab.h"
#include "clap.h"

int sim_bmark(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(rsiz,    int,     5,            "CA rule size");
	CLAP_CARG(fsiz,    int,     5,            "CA filter rule size");
	CLAP_CARG(n,       size_t,  10,           "number of words");
	CLAP_CARG(I,       size_t,  1000,         "number of iterations");
	CLAP_CARG(S,       size_t,  1000,         "number of samples");
	CLAP_CARG(maxrot,  int,     100,          "maximum rotation");
	CLAP_CARG(seed,    ulong,   0,            "random seed (0 for unpredictable)");
	puts("---------------------------------------------------------------------------------------\n");

	double ts,te;

	// pseudo-random number generation

	mt_t rng;
	mt_seed(&rng,seed);

	// CA rule tables

	word_t* const rtab = rt_alloc(rsiz);
	word_t* const ftab = rt_alloc(fsiz);
	rt_randomise(rsiz,rtab,0.5,&rng);
	rt_randomise(fsiz,ftab,0.5,&rng);

	// allocate CA storage

	const size_t N = I*n;
	const size_t T = S*N;

	word_t* const cas = mw_alloc(T);
	word_t* const uas = mw_alloc(T);
	word_t* const fas = mw_alloc(T);

	word_t** ca = malloc(S*sizeof(word_t*));
	TEST_ALLOC(ca);
	word_t** ua = malloc(S*sizeof(word_t*));
	TEST_ALLOC(ua);
	word_t** fa = malloc(S*sizeof(word_t*));
	TEST_ALLOC(fa);

	for (size_t k=0; k<S; ++k) ca[k] = cas+k*n;
	for (size_t k=0; k<S; ++k) ua[k] = uas+k*n;
	for (size_t k=0; k<S; ++k) fa[k] = fas+k*n;

	// rotations

	int* b = malloc(S*sizeof(int));
	TEST_ALLOC(b);
	for (int* p=b; p<b+S; ++p) *p = RANDI(int,maxrot,&rng);

	// intialise CA

	for (size_t k=0; k<S; ++k) mw_randomise(n,ca[k],&rng);

	// run CAs

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	for (size_t k=0; k<S; ++k) ca_run(I,n,ca[k],NULL,rsiz,rtab,0);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("CA run    time = %8.6f\n",te-ts);

	// rotate CAs

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	for (size_t k=0; k<S; ++k) ca_rotl(I,n,ua[k],ca[k],b[k]);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("CA rotate time = %8.6f\n",te-ts);

	// filter CAs

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	for (size_t k=0; k<S; ++k) ca_filter(I,n,fa[k],ua[k],fsiz,ftab);
	te = (double)clock()/(double)CLOCKS_PER_SEC;
	printf("CA filter time = %8.6f\n",te-ts);

	newline;

	// clean up

	free(b);
	free(fa);
	free(ua);
	free(ca);
	free(fas);
	free(uas);
	free(cas);
	free(ftab);
	free(rtab);

	return EXIT_SUCCESS;
}
