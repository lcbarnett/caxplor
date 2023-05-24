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

	printf("rtfile = %s\n",rtfile);

	FILE* const rtfs = fopen(rtfile,"r");
	if (rtfs == NULL) PEEXIT("failed to open saved rtids file '%s'",rtfile);

	rtl_t* rule  = NULL;

	char* line = NULL;
	size_t len = 0;
	int nlines = 0;
	int res;
	ssize_t ilen;
	while ((ilen = getline(&line,&len,rtfs)) != -1) {

		++nlines;
		printf("\nline %d: %s",nlines,line);

		line[ilen-1] = '\0'; // strip trailing newline

		const char* token = strtok(line," ");
		printf("\t(%s)\n",token);
		int rsiz;
		sscanf(token,"%d",&rsiz);
		printf("\tread %d\n",rsiz);
		if (rsiz < 1) {
			printf("\tbad length\n");
			continue;
		}

		token = strtok(NULL," ");
		printf("\t(%s)\n",token);
		printf("\tread [%s]\n",token);

		word_t* const rtab = rt_alloc(rsiz);
		res = rt_sread_id(rsiz,rtab,token);
		if (res == 1) {
			printf("\trtid is wrong length\n");
			free(rtab);
			continue;
		}
		if (res == 2) {
			printf("\trtid contains non-hex characters\n");
			free(rtab);
			continue;
		}
		rtl_t* erule = rtl_find(rule,rsiz,rtab);
		if (erule == NULL) {
			rule = rtl_add(rule,rsiz);
			printf("\tnew rule at %p\n",rule);
			rt_copy(rsiz,rule->tab,rtab);
		}
		else {
			rule = erule;
			printf("\trule exists at %p\n",rule);
		}
		free(rtab);

		token = strtok(NULL," ");
		printf("\t(%s)\n",token);
		printf("\tread [%s]\n",token);

		if (token == NULL) continue; // no filter id

		token = strtok(NULL," ");
		printf("\t(%s)\n",token);
		int fsiz;
		sscanf(token,"%d",&fsiz);
		printf("\tread %d\n",fsiz);
		if (rsiz < 1) {
			printf("\tbad length\n");
			continue;
		}


    }

	printf("\nnlines = %d\n\n",nlines);

	free(line);

	rtl_free(rule);

	if (fclose(rtfs) == -1) PEEXIT("failed to close saved rtids file '%s'",rtfile);

	return EXIT_SUCCESS;
}
