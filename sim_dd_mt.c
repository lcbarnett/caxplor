#include <pthread.h>

#include "clap.h"
#include "rtab.h"

#define ofnlen 1000

typedef struct {
	rtl_t* rule;
	rtl_t* filt;
	char   ofname[ofnlen+1];
} tfilt_t;

typedef struct {
	int tnum;
	int nfint;
	int emmax;
	int eiff;
	int tmmax;
	int tiff;
	int tlag;
	tfilt_t* tfilt;
} targ_t;

void* compfun(void* arg)
{
	const pthread_t tpid = pthread_self();
	const targ_t* const targ = (targ_t*)arg;
	const int tnum  = targ->tnum;
	const int nfint = targ->nfint;

	printf("thread %2d (%zu) : %d filters : STARTED\n",tnum,tpid,nfint);

	const int emmax = targ->emmax;
	const int eiff  = targ->eiff;
	const int tmmax = targ->tmmax;
	const int tiff  = targ->tiff;
	const int tlag  = targ->tlag;
	const int hlen  = (emmax > tmmax ? emmax : tmmax)+1;

	const tfilt_t* const tfilt = targ->tfilt;

	double Hr[hlen];
	double Hf[hlen];
	double DD[hlen];

	for (int i=0; i<nfint; ++i) {

		const int           rsize  = tfilt[i].rule->size;
		const int           fsize  = tfilt[i].filt->size;
		const word_t* const rtab   = tfilt[i].rule->tab;
		const word_t* const ftab   = tfilt[i].filt->tab;
		const char*   const ofname = tfilt[i].ofname;

		printf("\tthread %2d : filter %2d of %2d : rule id = ",tnum,i,nfint);
		rt_print_id(rsize,rtab);
		printf(", filter id = ");
		rt_print_id(fsize,ftab);
		putchar('\n'); fflush(stdout);

		const int rfsize = rsize > fsize ? rsize : fsize;

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

		FILE* const dfs = fopen(ofname,"w");
		if (dfs == NULL) PEEXIT("thread %2d : failed to open output file \"%s\"\n",tnum,ofname);
		for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,Hr[m],Hf[m],DD[m]);
		if (fclose(dfs) == -1) PEEXIT("thread %2d : failed to close output file\n",tnum);

		printf("\tthread %2d : filter %2d of %2d : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f : written \"%s\"\n",tnum,i,nfint,Hr[emmax],Hf[emmax],DD[tmmax],ofname);
	}

	printf("thread %2d (%zu) : %d filters : FINISHED\n",tnum,tpid,nfint);


	pthread_exit(NULL);
}

int sim_dd_mt(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(irtfile,  cstr,   "saved2.rt",   "input rtids file");
	CLAP_CARG(emmax,    int,     20,           "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,     int,     1,            "advance before entropy");
	CLAP_CARG(tmmax,    int,     14,           "maximum sequence length for DD calculation");
	CLAP_CARG(tiff,     int,     0,            "advance before DD calculation");
	CLAP_CARG(tlag,     int,     1,            "lag for DD calculation");
	CLAP_CARG(nthreads, int,     4,            "number of threads");
	CLAP_CARG(odir,     cstr,   "/tmp",        "output file directory");
	puts("---------------------------------------------------------------------------------------\n");

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
	int ntfilts = 0;
	for (size_t k=0;k<nrules;++k) ntfilts += (int)nfilts[k];
	printf("number of rules = %zu, filters =",nrules);
	for (size_t k=0;k<nrules;++k) printf(" %zu",nfilts[k]);
	printf(" (total = %d)\n\n",ntfilts);

	const int nfpert = ntfilts/nthreads+1;
	printf("threads = %d\nfilters = %d\nfilters per thread = %d (%d)\n\n",nthreads,ntfilts,nfpert,ntfilts%nthreads);

	// initialize and set thread joinable

	pthread_t threads[nthreads];
	pthread_attr_t attr;

	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	// thread-independent parameters

	targ_t targ[nthreads];
	for (int tnum=0; tnum<nthreads; ++tnum) {
		targ[tnum].emmax = emmax;
		targ[tnum].eiff  = eiff;
		targ[tnum].tmmax = tmmax;
		targ[tnum].tiff  = tiff;
		targ[tnum].tlag  = tlag;
		targ[tnum].tfilt = malloc((size_t)nfpert*sizeof(tfilt_t));
	}

	// loop through rules/filters

	int nfint = 0;
	int tnum = 0;
	for (; rule != NULL; rule = rule->next) {
		for (; rule->filt != NULL; rule->filt = rule->filt->next) {
			targ[tnum].tfilt[nfint].rule = rule;
			targ[tnum].tfilt[nfint].filt = rule->filt;
			snprintf(targ[tnum].tfilt[nfint].ofname,ofnlen,"%s/cadd_",odir);
			rt_sprint_id(rule->size,rule->tab,ofnlen,targ[tnum].tfilt[nfint].ofname);
			strncat(targ[tnum].tfilt[nfint].ofname,"_",1);
			rt_sprint_id(rule->filt->size,rule->filt->tab,ofnlen,targ[tnum].tfilt[nfint].ofname);
			strncat(targ[tnum].tfilt[nfint].ofname,".dat",4);
			if (++nfint == nfpert) {
				targ[tnum].tnum  = tnum;
				targ[tnum].nfint = nfint;
				const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targ[tnum]);
				if (tres) PEEXIT("unable to create thread %d",tnum+1)
				nfint = 0;
				++tnum;
			}
		}
	}
	if (nfint > 0) {
		targ[tnum].tnum  = tnum;
		targ[tnum].nfint = nfint;
		const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targ[tnum]);
		if (tres) PEEXIT("unable to create thread %d",tnum+1)
	}
//return 0;
	// free attribute and wait for the other threads to complete

	pthread_attr_destroy(&attr);
	for (int tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_join(threads[tnum],NULL);
		if (tres) PEEXIT("unable to join thread %d",tnum+1)
	}

	// clean up

	for (int tnum=0; tnum<nthreads; ++tnum) free(targ[tnum].tfilt);

	pthread_exit(NULL);

	free(nfilts);
	rtl_free(rule);

	return EXIT_SUCCESS;
}
