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

	return EXIT_SUCCESS;
}
