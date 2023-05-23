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
	CLAP_CARG(n,       size_t,  10,           "number of words");
	CLAP_CARG(B,       int,     5,            "CA rule breadth");
	CLAP_CARG(rlam,    double,  0.6,          "CA rule lambda");
	CLAP_CARG(rseed,   ulong,   0,            "CA rule random seed (or 0 for unpredictable)");
	CLAP_CARG(iseed,   ulong,   0,            "initialisation random seed (0 for unpredictable)");
	CLAP_CARG(S,       size_t,  100000,       "Sample size");
	puts("---------------------------------------------------------------------------------------\n");

	double ts,te;

	mt_t rrng, irng;
	mt_seed(&rrng,rseed);
	mt_seed(&irng,iseed);

	word_t* const rtab = rt_alloc(B);
	rt_randomise(B,rtab,rlam,&rrng);

	word_t* const w  = mw_alloc(n);
	word_t* const w1 = mw_alloc(n);
	word_t* const w2 = mw_alloc(n);
	mw_randomise(n,w,&irng);

	if (B > 6) {
		const size_t m = rt_nwords(B);
		word_t* const rt = mw_alloc(m);
		for (size_t r=0;r<POW2(B);++r) {const int k = (int)r; PUTBIT(rt[k/B],k%B,rtab[r]);}

		mw_filter(n,w1,w,B,rtab);
		mw_filtx (n,w2,w,B,rt  );

		mw_print(n,w1); putchar('\n');
		mw_print(n,w2); putchar('\n');

		ts = (double)clock()/(double)CLOCKS_PER_SEC;
		for (size_t k=0; k<S; ++k) mw_filter(n,w1,w,B,rtab);
		te = (double)clock()/(double)CLOCKS_PER_SEC;
		printf("\nOld method = %8.6f\n",te-ts);

		ts = (double)clock()/(double)CLOCKS_PER_SEC;
		for (size_t k=0; k<S; ++k) mw_filtx(n,w2,w,B,rt);
		te = (double)clock()/(double)CLOCKS_PER_SEC;
		printf("New method = %8.6f\n",te-ts);

		free(rt);
	}
	else {
		word_t rt;
		for (size_t r=0;r<POW2(B);++r) PUTBIT(rt,r,rtab[r]);

		mw_filter(n,w1,w,B,rtab);
		mw_filt  (n,w2,w,B,rt  );

		mw_print(n,w1); putchar('\n');
		mw_print(n,w2); putchar('\n');

		ts = (double)clock()/(double)CLOCKS_PER_SEC;
		for (size_t k=0; k<S; ++k) mw_filter(n,w1,w,B,rtab);
		te = (double)clock()/(double)CLOCKS_PER_SEC;
		printf("\nOld method = %8.6f\n",te-ts);

		ts = (double)clock()/(double)CLOCKS_PER_SEC;
		for (size_t k=0; k<S; ++k) mw_filt(n,w2,w,B,rt);
		te = (double)clock()/(double)CLOCKS_PER_SEC;
		printf("New method = %8.6f\n",te-ts);

	}

	printf("\n%s\n\n",mw_equal(n,w1,w2)?"good":"bad");
	free(w2);
	free(w1);
	free(w);
	free(rtab);

	return EXIT_SUCCESS;
}
