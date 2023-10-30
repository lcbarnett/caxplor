#include <pthread.h>

#include "clap.h"
#include "rtab.h"

#define odbuflen (200)
static char odir[odbuflen];

typedef struct {
	size_t tnum;
	rtl_t* rule;
	rtl_t* filt;
	int    emmax;
	int    eiff;
	int    tmmax;
	int    tiff;
	int    tlag;
} targ_t;

void* compfun(void* arg)
{
	const targ_t* const targ = (targ_t*)arg;

	const size_t tnum  = targ->tnum+1;
	const int rsize = targ->rule->size;
	const int fsize = targ->filt->size;
	const word_t* const rtab = targ->rule->tab;
	const word_t* const ftab = targ->filt->tab;
	const int emmax  = targ->emmax;
	const int eiff   = targ->eiff;
	const int tmmax  = targ->tmmax;
	const int tiff   = targ->tiff;
	const int tlag   = targ->tlag;
	const int hlen   = (emmax > tmmax ? emmax : tmmax)+1;
	const int rfsize = rsize > fsize ? rsize : fsize;

	printf("thread %2zu : rule id = ",tnum);
	rt_print_id(rsize,rtab);
	printf(", filter id = ");
	rt_print_id(fsize,ftab);
	putchar('\n'); fflush(stdout);

	double Hr[hlen];
	double Hf[hlen];
	double DD[hlen];

	for (int m=0; m<hlen; ++m) Hr[m] = NAN;
	for (int m=0; m<hlen; ++m) Hf[m] = NAN;
	for (int m=0; m<hlen; ++m) DD[m] = NAN;
	for (int m=rsize; m<=emmax; ++m) {
		Hr[m] = rt_entro(rsize,rtab,m,eiff)/(double)m;
	}
	for (int m=fsize; m<=emmax; ++m) {
		Hf[m] = rt_entro(fsize,ftab,m,eiff)/(double)m;
	}
	for (int m=rfsize; m<=tmmax; ++m) {
		DD[m] = rt_trent1(rsize,rtab,fsize,ftab,m,tiff,tlag)/(double)m;
	}

	const size_t ofnlen = 1000;
	char ofname[ofnlen+1];
	snprintf(ofname,ofnlen,"%s/cadd_",odir);
	rt_sprint_id(rsize,rtab,ofnlen,ofname);
	strncat(ofname,"_",1);
	rt_sprint_id(fsize,ftab,ofnlen,ofname);
	strncat(ofname,".dat",4);

	FILE* const dfs = fopen(ofname,"w");
	if (dfs == NULL) PEEXIT("thread %2zu : failed to open output file \"%s\"\n",tnum,ofname);
	for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,Hr[m],Hf[m],DD[m]);
	if (fclose(dfs) == -1) PEEXIT("thread %2zu : failed to close output file\n",tnum);

	printf("thread %2zu : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f : written \"%s\"\n",tnum,Hr[emmax],Hf[emmax],DD[tmmax],ofname);

	pthread_exit(NULL);
}

int sim_dd_mt(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(irtfile, cstr,   "saved.rt",    "input rtids file");
	CLAP_CARG(emmax,   int,     20,           "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,    int,     1,            "advance before entropy");
	CLAP_CARG(tmmax,   int,     14,           "maximum sequence length for DD calculation");
	CLAP_CARG(tiff,    int,     0,            "advance before DD calculation");
	CLAP_CARG(tlag,    int,     1,            "lag for DD calculation");
	CLAP_CARG(outdir,  cstr,   "/tmp",        "output file directory");
	puts("---------------------------------------------------------------------------------------\n");

	strncpy(odir,outdir,odbuflen);

	// Read in rule/filter rtids

	ASSERT(irtfile[0] != '\0',"Must supply an input rtid file");
	printf("Reading rules and filters from '%s' ...\n",irtfile);
	FILE* const irtfs = fopen(irtfile,"r");
	if (irtfs == NULL) PEEXIT("failed to open input rtids file '%s'",irtfile);
	rtl_t* rule = rtl_fread(irtfs);
	ASSERT(rule != NULL,"No valid rtids found in input file!");
	if (fclose(irtfs) == -1) PEEXIT("failed to close input rtids file '%s'",irtfile);

	// Calculate number of rules/filters (number of threads = total number of filters)

	size_t nrules;
	size_t* const nfilts = rtl_nitems(rule,&nrules);
	size_t ntfilts = 0;
	for (size_t k=0;k<nrules;++k) ntfilts += nfilts[k];
	printf("number of rules = %zu, filters =",nrules);
	for (size_t k=0;k<nrules;++k) printf(" %zu",nfilts[k]);
	printf(" (total = %zu)\n\n",ntfilts);

	// initialize and set thread joinable

	pthread_t threads[ntfilts];
	pthread_attr_t attr;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	// set up parameter struct array

	targ_t targ[ntfilts];
	for (size_t tnum=0; tnum<ntfilts; ++tnum) {
		targ[tnum].tnum  = tnum;
		targ[tnum].emmax = emmax;
		targ[tnum].eiff  = eiff;
		targ[tnum].tmmax = tmmax;
		targ[tnum].tiff  = tiff;
		targ[tnum].tlag  = tlag;
	}

	// kick off computation threads
	size_t tnum = 0;
	for (; rule != NULL; rule = rule->next) {
		for (; rule->filt != NULL; rule->filt = rule->filt->next,++tnum) {
			targ[tnum].rule  = rule;
			targ[tnum].filt  = rule->filt;
			const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targ[tnum]);
			if (tres) {
				fprintf(stderr,"Error: unable to create thread %zu",tnum);
				exit(EXIT_FAILURE);
			}
		}
	}

	// free attribute and wait for the other threads
	pthread_attr_destroy(&attr);
	for (size_t tnum=0; tnum<ntfilts; ++tnum) {
		const int tres = pthread_join(threads[tnum],NULL);
		if (tres) {
			fprintf(stderr,"Error: unable to join thread %zu",tnum);
			exit(EXIT_FAILURE);
		}
	}

	pthread_exit(NULL);

	free(nfilts);
	rtl_free(rule);

	return EXIT_SUCCESS;
}
