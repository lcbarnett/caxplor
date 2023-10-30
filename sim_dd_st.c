#include "clap.h"
#include "rtab.h"

int sim_dd_st(int argc, char* argv[])
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
	CLAP_CARG(odir,    cstr,   "/tmp",        "output file directory");
	puts("---------------------------------------------------------------------------------------\n");

	// Read in rule/filter rtids

	ASSERT(irtfile[0] != '\0',"Must supply an input rtid file");
	printf("Reading rules and filters from '%s' ...\n",irtfile);
	FILE* const irtfs = fopen(irtfile,"r");
	if (irtfs == NULL) PEEXIT("failed to open input rtids file '%s'",irtfile);
	rtl_t* rule = rtl_fread(irtfs);
	ASSERT(rule != NULL,"No valid rtids found in input file!");
	if (fclose(irtfs) == -1) PEEXIT("failed to close input rtids file '%s'",irtfile);

	size_t nrules;
	size_t* const nfilts = rtl_nitems(rule,&nrules);
	printf("number of rules = %zu, filters =",nrules);
	for (size_t k=0;k<nrules;++k) printf(" %zu",nfilts[k]);
	puts("\n");


	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	double Hr[hlen];
	double Hf[hlen];
	double DD[hlen];

	const size_t ofnlen = 999;
	char ofname[ofnlen+1];

	// Run through rule/filter rtids lists calculating entropy and DD
	for (; rule != NULL; rule = rule->next) {
		printf("rule id = "); rt_print_id(rule->size,rule->tab);
		for (int m=0; m<hlen; ++m) Hr[m] = NAN;
		for (int m=rule->size; m<=emmax; ++m) {
			Hr[m] = rt_entro(rule->size,rule->tab,m,eiff)/(double)m;
			putchar('.'); fflush(stdout);
		}
		printf(" rule entropy ≈ %8.6f\n",Hr[emmax]);
		for (; rule->filt != NULL; rule->filt = rule->filt->next) {
			printf("\tfilter id = "); rt_print_id(rule->filt->size,rule->filt->tab);
			for (int m=0; m<hlen; ++m) Hf[m] = NAN;
			for (int m=rule->filt->size; m<=emmax; ++m) {
				Hf[m] = rt_entro(rule->filt->size,rule->filt->tab,m,eiff)/(double)m;
				putchar('.'); fflush(stdout);
			}
			printf(" filter entropy ≈ %8.6f",Hf[emmax]); fflush(stdout);
			const int mmin = rule->size > rule->filt->size ? rule->size : rule->filt->size;
			for (int m=0; m<hlen; ++m) DD[m] = NAN;
			for (int m=mmin; m<=tmmax; ++m) {
				DD[m] = rt_trent1(rule->size,rule->tab,rule->filt->size,rule->filt->tab,m,tiff,tlag)/(double)m;
				putchar('.'); fflush(stdout);
			}
			printf(" DD ≈ %8.6f",DD[tmmax]);

			snprintf(ofname,ofnlen,"%s/cadd_",odir);
			rt_sprint_id(rule->size,rule->tab,ofnlen,ofname);
			strncat(ofname,"_",1);
			rt_sprint_id(rule->filt->size,rule->filt->tab,ofnlen,ofname);
			strncat(ofname,".dat",4);

			FILE* const dfs = fopen(ofname,"w");
			if (dfs == NULL) PEEXIT("failed to open output file \"%s\"\n",ofname);
			for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,Hr[m],Hf[m],DD[m]);
			if (fclose(dfs) == -1) PEEXIT("failed to close output file\n");
			printf(" : written \"%s\"\n",ofname);
		}
	}

	free(nfilts);
	rtl_free(rule);

	return EXIT_SUCCESS;
}
