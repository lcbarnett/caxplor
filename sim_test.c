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
	CLAP_CARG(n,       int,     42,           "a variable");
	puts("---------------------------------------------------------------------------------------\n");

	printf("n = %d\n",n);

	return EXIT_SUCCESS;
}
