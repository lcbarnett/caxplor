#include <time.h>

#include "ca.h"
#include "clap.h"
#include "utils.h"

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(rsiz,    int,     5,            "CA rule size");
	CLAP_VARG(fsiz,    int,     0,            "filter rule size (or 0 for same as rule size)");
	CLAP_CARG(rtfile,  cstr,   "saved.rt",    "saved rtids file name");
	puts("---------------------------------------------------------------------------------------\n");

	fsiz = fsiz == 0 ? rsiz : fsiz;

	printf("rsiz = %d, fsiz = %d, rtfile = %s\n\n",rsiz,fsiz,rtfile);

	FILE* const rtfs = fopen(rtfile,"r");
	if (rtfs == NULL) PEEXIT("failed to open saved rtids file '%s'",rtfile);

	char* line = NULL;
	size_t len = 0;
	int nlines = 0;
	while (getline(&line,&len,rtfs) != -1) {
		++nlines;
		printf("line %d: %s",nlines,line);

		const char* token = strtok(line," ");
		printf("\t(%s)\n",token);
		int BB;
		sscanf(token,"%d",&BB);
		printf("\tread %d\n",BB);

		token = strtok(NULL," ");
		printf("\t(%s)\n",token);
		printf("\tread [%s]\n",token);

    }

	printf("\nnlines = %d\n\n",nlines);

	free(line);

	if (fclose(rtfs) == -1) PEEXIT("failed to close saved rtids file '%s'",rtfile);

	return EXIT_SUCCESS;
}
