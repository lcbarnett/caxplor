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
	CLAP_CARG(rtfile,  cstr,   "saved.rt",    "saved rtids file name");
	puts("---------------------------------------------------------------------------------------\n");

	FILE* const rtfs = fopen(rtfile,"r");
	if (rtfs == NULL) PEEXIT("failed to open saved rtids file '%s'",rtfile);

	rtl_t* rule = rtl_fread(rtfs);

	rtl_free(rule);

	if (fclose(rtfs) == -1) PEEXIT("failed to close saved rtids file '%s'",rtfile);

	return EXIT_SUCCESS;
}
