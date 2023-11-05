#include <pthread.h>

#include "clap.h"
#include "rtab.h"

typedef struct {
	rtl_t*  rule;
	rtl_t*  filt;
	double* Hr;
	double* Hf;
	double* DD;
} tfarg_t;

typedef struct {
	int tnum;
	int nfint;
	int emmax;
	int eiff;
	int tmmax;
	int tiff;
	int tlag;
	tfarg_t* tfargs;
} targ_t;

void* compfun(void* arg);

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

	targ_t targs[nthreads];
	for (int tnum=0; tnum<nthreads; ++tnum) {
		targs[tnum].emmax  = emmax;
		targs[tnum].eiff   = eiff;
		targs[tnum].tmmax  = tmmax;
		targs[tnum].tiff   = tiff;
		targs[tnum].tlag   = tlag;
		targs[tnum].tfargs = malloc((size_t)nfpert*sizeof(tfarg_t));
		TEST_ALLOC(targs[tnum].tfargs);
	}

	// loop through rules/filters, setting up thread-dependent parameters

	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	int tnum  = 0;
	int nfint = 0;
	for (rtl_t* r = rtl_init(rule); r != NULL; r = r->next) {
		for (rtl_t* f = r->filt; f != NULL; f = f->next) {
			targ_t*  const targ  = &targs[tnum];
			tfarg_t* const tfarg = &targ->tfargs[nfint];
			tfarg->rule = r;
			tfarg->filt = f;
			tfarg->Hr = malloc((size_t)hlen*sizeof(double));
			tfarg->Hf = malloc((size_t)hlen*sizeof(double));
			tfarg->DD = malloc((size_t)hlen*sizeof(double));
			if (++nfint == nfpert) {
				targ->tnum  = tnum;
				targ->nfint = nfint;
				++tnum;
				nfint = 0;
			}
		}
	}
	if (nfint > 0) {
		targ_t* const targ = &targs[tnum];
		targ->tnum  = tnum;
		targ->nfint = nfint;
	}

	// set up computational threads as joinable

	pthread_t threads[nthreads];
	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setdetachstate(&attr,PTHREAD_CREATE_JOINABLE);

	// create threads

	for (int tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_create(&threads[tnum],&attr,compfun,(void*)&targs[tnum]);
		PASSERT(tres == 0,"unable to create thread %d",tnum+1)
	}

	// free attribute and wait for computational threads to complete

	pthread_attr_destroy(&attr);
	for (int tnum=0; tnum<nthreads; ++tnum) {
		const int tres = pthread_join(threads[tnum],NULL);
		PASSERT(tres == 0,"unable to join thread %d",tnum+1)
	}

	// write out results

	const size_t ofnlen = strlen(odir)+10;
	char ofname[ofnlen];
	snprintf(ofname,ofnlen,"%s/cadd.dat",odir);
	printf("\nWriting results to \"%s\"... ",ofname);
	fflush(stdout);
	FILE* const dfs = fopen(ofname,"w");
	PASSERT(dfs != NULL,"Failed to open output file \"%s\"\n",ofname);
	for (int tnum=0; tnum<nthreads; ++tnum) {
		targ_t* const targ = &targs[tnum];
		for (int i=0; i<targ->nfint; ++i) {
			tfarg_t* const tfarg = &targ->tfargs[i];
			fprintf(dfs,"# rule id = ");
			rt_fprint_id(tfarg->rule->size,tfarg->rule->tab,dfs);
			fprintf(dfs,", filter id = ");
			rt_fprint_id(tfarg->filt->size,tfarg->filt->tab,dfs);
			fputc('\n',dfs);
			for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,tfarg->Hr[m],tfarg->Hf[m],tfarg->DD[m]);
			fputs("\n",dfs);
		}
	}
	if (fclose(dfs) == -1) PEEXIT("Failed to close output file \"%s\"\n",ofname);
	puts("done");

	// clean up

	for (int tnum=0; tnum<nthreads; ++tnum) {
		targ_t* const targ = &targs[tnum];
		for (int i=0; i<targ->nfint; ++i) {
			tfarg_t* const tfarg = &targ->tfargs[i];
			free(tfarg->DD);
			free(tfarg->Hf);
			free(tfarg->Hr);
		}
		free(targ->tfargs);
	}

	rtl_free(rule);

	return EXIT_SUCCESS;
}

void* compfun(void* arg)
{
	const pthread_t tpid = pthread_self();
	const targ_t* const targs = (targ_t*)arg;
	const int tnum  = targs->tnum;
	const int nfint = targs->nfint;

	printf("thread %2d (%zu) : %d filters : STARTED\n",tnum+1,tpid,nfint);
	fflush(stdout);

	const int emmax = targs->emmax;
	const int eiff  = targs->eiff;
	const int tmmax = targs->tmmax;
	const int tiff  = targs->tiff;
	const int tlag  = targs->tlag;
	const int hlen  = (emmax > tmmax ? emmax : tmmax)+1;

	const tfarg_t* const tfargs = targs->tfargs;

	for (int i=0; i<nfint; ++i) {

		const int           rsize = tfargs[i].rule->size;
		const int           fsize = tfargs[i].filt->size;
		const word_t* const rtab  = tfargs[i].rule->tab;
		const word_t* const ftab  = tfargs[i].filt->tab;
		double*       const Hr    = tfargs[i].Hr;
		double*       const Hf    = tfargs[i].Hf;
		double*       const DD    = tfargs[i].DD;

		const int rfsize = rsize > fsize ? rsize : fsize;

		for (int m=0; m<hlen; ++m) Hr[m] = NAN;
		for (int m=0; m<hlen; ++m) Hf[m] = NAN;
		for (int m=0; m<hlen; ++m) DD[m] = NAN;
		for (int m=rsize;  m<=emmax; ++m) Hr[m] = rt_entro(rsize,rtab,m,eiff)/(double)m;
		for (int m=fsize;  m<=emmax; ++m) Hf[m] = rt_entro(fsize,ftab,m,eiff)/(double)m;
		for (int m=rfsize; m<=tmmax; ++m) DD[m] = rt_trent1(rsize,rtab,fsize,ftab,m,tiff,tlag)/(double)m;

		flockfile(stdout); // prevent another thread butting in!
		printf("\tthread %2d : filter %2d of %2d : rule id = ",tnum+1,i+1,nfint);
		rt_print_id(rsize,rtab);
		printf(", filter id = ");
		rt_print_id(fsize,ftab);
		printf(" : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f\n",Hr[emmax],Hf[emmax],DD[tmmax]);
		fflush(stdout);
		funlockfile(stdout);
	}

	printf("thread %2d : FINISHED\n",tnum+1);
	fflush(stdout);

	pthread_exit(NULL);
}
