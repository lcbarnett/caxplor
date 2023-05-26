#include <time.h>

#include "ca.h"
#include "rtab.h"
#include "clap.h"
#include "utils.h"


rtl_t* rtl_read_file(FILE* rtfs)
{

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
		int rsiz;
		sscanf(token,"%d",&rsiz);
		printf("\trsiz: (%s) read %d\n",token,rsiz);
		if (rsiz < 1) {
			printf("\tbad length\n");
			continue;
		}

		token = strtok(NULL," ");
		printf("\trtid: read [%s]\n",token);

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
		rtl_t* rrule = rtl_find(rule,rsiz,rtab);
		if (rrule == NULL) {
			rule = rtl_add(rule,rsiz);
			printf("\tnew rule at %p\n",rule);
			rt_copy(rsiz,rule->tab,rtab);
		}
		else {
			rule = rrule;
			printf("\trule exists at %p\n",rule);
		}
		free(rtab);

		token = strtok(NULL," ");
		if (token == NULL) { // no filter id
			printf("\t(%s) no filter id\n",token);
			continue;
		}

		int fsiz;
		sscanf(token,"%d",&fsiz);
		printf("\tfsiz: (%s) read %d\n",token,fsiz);
		if (fsiz < 1) {
			printf("\tbad length\n");
			continue;
		}

		token = strtok(NULL," ");
		printf("\tftid: read [%s]\n",token);

		word_t* const ftab = rt_alloc(fsiz);
		res = rt_sread_id(fsiz,ftab,token);
		if (res == 1) {
			printf("\tftid is wrong length\n");
			free(ftab);
			continue;
		}
		if (res == 2) {
			printf("\tftid contains non-hex characters\n");
			free(ftab);
			continue;
		}
		rtl_t* frule = rtl_find(rule->filt,fsiz,ftab);
		if (frule == NULL) {
			rule->filt = rtl_add(rule->filt,fsiz);
			printf("\tnew filter rule at %p\n",rule);
			rt_copy(fsiz,rule->filt->tab,ftab);
		}
		else {
			rule->filt = frule;
			printf("\tfilter rule exists at %p\n",rule->filt);
		}
		free(ftab);

		token = strtok(NULL," ");
		if (token != NULL) { // extra token
			printf("\t(%s) extra tokens!\n",token);
			continue;
		}
    }

	printf("\nnlines = %d\n\n",nlines);

	free(line);

	while (rule->prev != NULL) rule = rule->prev; // rewind to beginning of list

	return rule;
}

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

	rtl_t* rule = rtl_read_file(rtfs);

	rtl_free(rule);

	if (fclose(rtfs) == -1) PEEXIT("failed to close saved rtids file '%s'",rtfile);

	return EXIT_SUCCESS;
}
