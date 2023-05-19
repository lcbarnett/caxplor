#include <time.h>

#include "ca.h"
#include "clap.h"

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(B,       int,     5,            "CA rule breadth");
	CLAP_CARG(rbias,   double,  0.5,          "CA rule bias");
	CLAP_CARG(b,       size_t,  7,            "number of bits to set");
	CLAP_CARG(n,       size_t,  9,            "number of words");
	CLAP_CARG(S,       size_t,  1000,         "number of samples");
	CLAP_CARG(maxrot,  int,     100,          "maximum rotation");
	CLAP_CARG(rseed,   ulong,   0,            "random seed (0 for unpredictable)");
	CLAP_CARG(wseed,   ulong,   0,            "random seed (0 for unpredictable)");
	puts("---------------------------------------------------------------------------------------\n");

	mt_t rng;
	mt_seed(&rng,rseed);

	word_t* const rtab = rt_alloc(B);
	rt_randomb(B,rtab,b,&rng);

	rt_print(B,rtab);

	free(rtab);


/*
struct timespec res;
int ret;
ret = clock_getres(CLOCK_REALTIME,&res);
PASSERT(ret == 0,"clock_getres failed\n");
printf("CLOCK_REALTIME : %ld\n",res.tv_nsec);
ret = clock_getres(CLOCK_REALTIME_COARSE,&res);
PASSERT(ret == 0,"clock_getres failed\n");
printf("CLOCK_REALTIME_COARSE : %ld\n",res.tv_nsec);
ret = clock_getres(CLOCK_MONOTONIC,&res);
PASSERT(ret == 0,"clock_getres failed\n");
printf("CLOCK_MONOTONIC : %ld\n",res.tv_nsec);
ret = clock_getres(CLOCK_PROCESS_CPUTIME_ID,&res);
PASSERT(ret == 0,"clock_getres failed\n");
printf("CLOCK_PROCESS_CPUTIME_ID : %ld\n",res.tv_nsec);
ret = clock_getres(CLOCK_THREAD_CPUTIME_ID,&res);
PASSERT(ret == 0,"clock_getres failed\n");
printf("CLOCK_THREAD_CPUTIME_ID : %ld\n",res.tv_nsec);
return 0;
*/

/*
	mt_t rrng, wrng;
	mt_seed(&rrng,rseed);
	mt_seed(&wrng,wseed);

	word_t* const rtab = rt_alloc(B);
	rt_randomise(B,rtab,rbias,&rrng);

	size_t Sn = S*n;

	word_t*   words = malloc(Sn*sizeof(word_t));
	word_t**  word  = malloc(S*sizeof(word_t*));
	for (size_t s=0; s<S; ++s) word[s] = words+s*n;

	word_t*   fwords = malloc(Sn*sizeof(word_t));
	word_t**  fword  = malloc(S*sizeof(word_t*));
	for (size_t s=0; s<S; ++s) fword[s] = fwords+s*n;

	word_t*   f1words = malloc(Sn*sizeof(word_t));
	word_t**  f1word  = malloc(S*sizeof(word_t*));
	for (size_t s=0; s<S; ++s) f1word[s] = f1words+s*n;

	for (word_t* p=words; p<words+Sn; ++p) *p = mt_uint(&wrng);

	double ts,ten,teo;

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	for (size_t s=0; s<S; ++s) mw_filter_old(n,f1word[s],word[s],B,rtab);
	teo = (double)clock()/(double)CLOCKS_PER_SEC-ts;
	printf("old filter = %8.6f\n",teo);

	ts = (double)clock()/(double)CLOCKS_PER_SEC;
	for (size_t s=0; s<S; ++s) mw_filter(n,fword[s],word[s],B,rtab);
	ten = (double)clock()/(double)CLOCKS_PER_SEC-ts;
	printf("new filter = %8.6f\n",ten);
	newline;
	printf("change = %4.2f%%\n",100.0*((teo-ten)/teo));
	newline;

	size_t e = 0;
	for (size_t s=0; s<S; ++s) {
		if (!mw_equal(n,fword[s],f1word[s])) ++e;
	}
	printf("errors = %lu\n",e);
	newline;

	free(f1words);
	free(fwords);
	free(words);
	free(rtab);
*/
	return EXIT_SUCCESS;
}
