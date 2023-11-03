#include <pthread.h>

#include "clap.h"
#include "rtab.h"

#define ofnmaxlen 1000

typedef struct {
	rtl_t* rule;
	rtl_t* filt;
	char*  ofname;
} tfarg_t;

typedef struct {
	int tnum;
	int nfint;
	int emmax;
	int eiff;
	int tmmax;
	int tiff;
	int tlag;
	tfarg_t* tfarg;
} targ_t;

void* compfun(void* arg)
{
	const pthread_t tpid = pthread_self();
	const targ_t* const targ = (targ_t*)arg;
	const int tnum  = targ->tnum;
	const int nfint = targ->nfint;

	printf("thread %2d (%zu) : %d filters : STARTED\n",tnum+1,tpid,nfint);
	fflush(stdout);

	const int emmax = targ->emmax;
	const int eiff  = targ->eiff;
	const int tmmax = targ->tmmax;
	const int tiff  = targ->tiff;
	const int tlag  = targ->tlag;
	const int hlen  = (emmax > tmmax ? emmax : tmmax)+1;

	const tfarg_t* const tfarg = targ->tfarg;

	double Hr[hlen];
	double Hf[hlen];
	double DD[hlen];

	for (int i=0; i<nfint; ++i) {

		const int           rsize  = tfarg[i].rule->size;
		const int           fsize  = tfarg[i].filt->size;
		const word_t* const rtab   = tfarg[i].rule->tab;
		const word_t* const ftab   = tfarg[i].filt->tab;
		const char*   const ofname = tfarg[i].ofname;

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
		PASSERT(dfs != NULL,"thread %2d : failed to open output file \"%s\"\n",tnum,ofname);
		for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,Hr[m],Hf[m],DD[m]);
		if (fclose(dfs) == -1) PEEXIT("thread %2d : failed to close output file\n",tnum);

		flockfile(stdout); // prevent another thread butting in!
		printf("\tthread %2d : filter %2d of %2d : rule id = ",tnum+1,i+1,nfint);
		rt_print_id(rsize,rtab);
		printf(", filter id = ");
		rt_print_id(fsize,ftab);
		printf(" : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f : written \"%s\"\n",Hr[emmax],Hf[emmax],DD[tmmax],ofname);
		fflush(stdout);
		funlockfile(stdout);
	}

	printf("thread %2d : FINISHED\n",tnum+1);
	fflush(stdout);

	pthread_exit(NULL);
}

size_t set_ofname(char* const ofname, const rtl_t* const r, const rtl_t* const f, const char* const odir)
{
	snprintf(ofname,ofnmaxlen,"%s/cadd_",odir);
	rt_sprint_id(r->size,r->tab,ofnmaxlen,ofname);
	strncat(ofname,"_",1);
	rt_sprint_id(f->size,f->tab,ofnmaxlen,ofname);
	strncat(ofname,".dat",4);
	return strlen(ofname);
}

int sim_dd_mt(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name      type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(irtfile,  cstr,   "saved.rt",    "input rtids file");
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
	PASSERT(irtfs != NULL,"failed to open input rtids file");
	rtl_t* rule = rtl_fread(irtfs);
	ASSERT(rule != NULL,"No valid rtids found in input file!");
	const int fres = fclose(irtfs);
	PASSERT(fres == 0,"failed to close input rtids file");

	// Calculate number of rules/filters

	int nrules, nfilts;
	int* const nfperr = rtl_nitems(rule,&nrules,&nfilts);
	printf("\nrules   = %d\n",nrules);
	printf("filters = %d :",nfilts);
	for (int k=0;k<nrules;++k) printf(" %d",nfperr[k]);
	putchar('\n');
	const int nfpert  = nfilts/nthreads + (nfilts%nthreads ? 1 : 0);
	const int nfpertl = nfilts-nfpert*(nthreads-1);
	printf("threads = %d\nfilters per thread = %d (last = %d)\n\n",nthreads,nfpert,nfpertl);
	fflush(stdout);
	ASSERT(nfpertl > 0,"Too many threads!?");
	free(nfperr);

	// thread-independent parameters

	targ_t targ[nthreads];
	for (int tnum=0; tnum<nthreads; ++tnum) {
		targ[tnum].emmax = emmax;
		targ[tnum].eiff  = eiff;
		targ[tnum].tmmax = tmmax;
		targ[tnum].tiff  = tiff;
		targ[tnum].tlag  = tlag;
		targ[tnum].tfarg = malloc((size_t)nfpert*sizeof(tfarg_t));
		TEST_ALLOC(targ[tnum].tfarg);
	}

	// loop through rules/filters, setting up thread-dependent parameters

	char ofname[ofnmaxlen+1];

	int tnum  = 0;
	int nfint = 0;
	for (rtl_t* r = rtl_init(rule); r != NULL; r = r->next) {
		for (rtl_t* f = r->filt; f != NULL; f = f->next) {
			tfarg_t* const tfarg =  &targ[tnum].tfarg[nfint];
			tfarg->rule = r;
			tfarg->filt = f;
			const size_t ofnlen = set_ofname(ofname,r,f,odir);
			tfarg->ofname = malloc(ofnlen+1);
			TEST_ALLOC(tfarg->ofname);
			strncpy(tfarg->ofname,ofname,ofnlen+1);
			if (++nfint == nfpert) {
				targ[tnum].tnum  = tnum;
				targ[tnum].nfint = nfint;
				++tnum;
				nfint = 0;
			}
		}
	}
	if (nfint > 0) {
		targ[tnum].tnum  = tnum;
		targ[tnum].nfint = nfint;
	}

	// set up computational threads as joinable

	pthread_t threads[nthreads];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	// create threads

	for (int tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targ[tnum]);
		PASSERT(tres == 0,"unable to create thread %d",tnum+1)
	}

	// free attribute and wait for computational threads to complete

	pthread_attr_destroy(&attr);
	for (int tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_join(threads[tnum],NULL);
		PASSERT(tres == 0,"unable to join thread %d",tnum+1)
	}

	// clean up

	for (int tnum=0; tnum<nthreads; ++tnum) {
		for (int i=0; i<targ[tnum].nfint; ++i) free(targ[tnum].tfarg[i].ofname);
		free(targ[tnum].tfarg);
	}

	rtl_free(rule);

	return EXIT_SUCCESS;
}
