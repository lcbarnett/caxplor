#include "clap.h"
#include "rtab.h"

int sim_dd(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(irtfile, cstr,   "",            "input rtids file (empty to start with random rtid)");
	CLAP_CARG(emmax,   int,     20,           "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,    int,     1,            "advance before entropy");
	CLAP_CARG(tmmax,   int,     14,           "maximum sequence length for DD calculation");
	CLAP_CARG(tiff,    int,     0,            "advance before DD calculation");
	CLAP_CARG(tlag,    int,     1,            "lag for DD calculation");
	puts("---------------------------------------------------------------------------------------\n");

	// Read in rule/filter rtids

	ASSERT(irtfile[0] != '\0',"Must supply an input rtid file");
	printf("Reading rules and filters from '%s' ...\n",irtfile);
	FILE* const irtfs = fopen(irtfile,"r");
	if (irtfs == NULL) PEEXIT("failed to open input rtids file '%s'",irtfile);
	rtl_t* rule = rtl_fread(irtfs);
	ASSERT(rule != NULL,"No valid rtids found in input file!");
	if (fclose(irtfs) == -1) PEEXIT("failed to close input rtids file '%s'",irtfile);
	printf("Done\n\n");

	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	double H [hlen];
	double Hf[hlen];
	double Tf[hlen];

	// Run through rule/filter rtids calculating te
	for (; rule != NULL; rule = rule->next) {
		for (int m=0; m<hlen; ++m) H[m] = NAN;
		for (int m=rule->size; m<=emmax; ++m) H[m]  = rt_entro(rule->size,rule->tab,m,eiff)/(double)m;
		for (; rule->filt != NULL; rule->filt = rule->filt->next) {
			printf("rule id = "); rt_print_id(rule->size,rule->tab);
			printf(", filter id = "); rt_print_id(rule->filt->size,rule->filt->tab);
			for (int m=0; m<hlen; ++m) Hf[m] = NAN;
			for (int m=rule->filt->size; m<=emmax; ++m) Hf[m]  = rt_entro(rule->filt->size,rule->filt->tab,m,eiff)/(double)m;
			const int mmin = rule->size > rule->filt->size ? rule->size : rule->filt->size;
			for (int m=0; m<hlen; ++m) Tf[m] = NAN;
			for (int m=mmin; m<=tmmax; ++m) Tf[m] = rt_trent1(rule->size,rule->tab,rule->filt->size,rule->filt->tab,m,tiff,tlag)/(double)m;
			printf(" rule entropy = %8.6f, filter entropy = %8.6f, DD = %8.6f\n",H[emmax],Hf[emmax],Tf[tmmax]);
		}
	}

	return EXIT_SUCCESS;
}
