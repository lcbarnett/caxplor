#include <time.h>

#include "ca.h"
#include "rtab.h"
#include "clap.h"
#include "utils.h"

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(len,     size_t,  5,            "number of hex chars");
	puts("---------------------------------------------------------------------------------------\n");

	for (size_t len = 0; len <= 128; ++len) {
		const int size = rt_hexsize(len);
		printf("len = %3lu, size = % d\n",len,size);
	}

	return EXIT_SUCCESS;
}
