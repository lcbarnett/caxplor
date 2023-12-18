#include <pthread.h>
#include <stdio.h>

#include "clap.h"
#include "rtab.h"

typedef struct {
	word_t* rtab;
	word_t* ftab;
	double* Hr;
	double* Hf;
	double* DD;
} tfarg_t;

typedef struct {
	size_t    tnum;
	size_t    nfpert;
	int       rsize;
	int       fsize;
	int       emmax;
	int       eiff;
	int       tmmax;
	int       tiff;
	int       tlag;
	uint64_t* ebuf;
	uint64_t* tbuf;
	tfarg_t*  tfargs;
} targ_t;

static void* compfun(void* arg);

int sim_ddr(int argc, char* argv[], int info)
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name      type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(rsize,    int,     5,             "CA rule size");
	CLAP_CARG(rlam,     double,  0.6,           "CA rule lambda");
	CLAP_CARG(rseed,    ulong,   0,             "CA rule random seed (or 0 for unpredictable)");
	CLAP_CARG(fsize,    int,     5,             "filter rule size");
	CLAP_CARG(flammin,  double,  0.6,           "filter rule lambda range minimum");
	CLAP_CARG(flammax,  double,  0.9,           "filter rule lambda range maximum");
	CLAP_CARG(flamres,  size_t,  10,            "filter rule lambda range resolution");
	CLAP_CARG(fseed,    ulong,   0,             "filter rule random seed (0 for unpredictable)");
	CLAP_CARG(emmax,    int,     20,            "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,     int,     1,             "advance before entropy");
	CLAP_CARG(tmmax,    int,     14,            "maximum sequence length for DD calculation");
	CLAP_CARG(tiff,     int,     0,             "advance before DD calculation");
	CLAP_CARG(tlag,     int,     1,             "lag for DD calculation");
	CLAP_CARG(nthreads, size_t,  4,             "number of threads");
	CLAP_CARG(nfpert,   size_t,  10,            "number of rules/filters per thread");
	CLAP_CARG(odir,     cstr,   "/tmp",         "output file directory");
	CLAP_CARG(jobidx,   cstr,   "LSB_JOBINDEX", "job index");
	puts("---------------------------------------------------------------------------------------\n");

	// set filter lambda parameter from job index

	size_t jnum;
	const char* const jobidxenv = getenv(jobidx);
	if (jobidxenv == NULL) {
		jnum = (size_t)atoi(jobidx);
		printf("*** Job number (from command line) %zu : ",jnum);
	}
	else {
		jnum = (size_t)atoi(jobidxenv);
		printf("*** Job number (from environmental variable %s) %zu :",jobidx,jnum);
	}
	fflush(stdout);
	ASSERT(jnum > 0 && jnum <= flamres, "Bad job index (%zu)",jnum);
	const double flam = flammin+(double)(jnum-1)*((flammax-flammin)/((double)(flamres-1)));
	printf("filter lambda = %g\n\n",flam);

	// buffer sizes for heap memory allocation

	const size_t rlen  = POW2(rsize);
	const size_t flen  = POW2(fsize);
	const size_t hlen  = (size_t)(emmax > tmmax ? emmax : tmmax)+1;
	const size_t eblen = POW2(emmax);
	const size_t tblen = POW2(2*tmmax);

	const unsigned long minmem =
		nthreads*nfpert*(rlen+flen)*sizeof(word_t) +
		nthreads*nfpert*3*hlen*sizeof(double) +
		nthreads*nfpert*sizeof(tfarg_t) +
		nthreads*(eblen+tblen)*sizeof(uint64_t);

	TEST_RAM(minmem);

	const double mmMb = (double)minmem/1000.0/1000.0;
	const double mmGb = mmMb/1000.0;
	printf("*** Dynamic memory > %.0fMb = %.2fGb\n\n",mmMb,mmGb);

	if (info) return EXIT_SUCCESS; // display some info and return

	// pseudo-random number generators

	mt_t rrng, frng;
	mt_seed(&rrng,rseed);
	mt_seed(&frng,fseed);

	// allocate buffers for random rules and filters

	printf("*** Allocating storage\n\n");

	word_t* const  rbuf = malloc(nthreads*nfpert*rlen*sizeof(word_t));
	TEST_ALLOC(rbuf);

	word_t* const  fbuf = malloc(nthreads*nfpert*flen*sizeof(word_t));
	TEST_ALLOC(fbuf);

	// allocate work buffers for entropy and DD computation

	uint64_t* const ebuf = malloc(nthreads*eblen*sizeof(uint64_t));
	TEST_ALLOC(ebuf);

	uint64_t* const tbuf = malloc(nthreads*tblen*sizeof(uint64_t));
	TEST_ALLOC(tbuf);

	// allocate storage buffers for entropy and DD results

	double* const  Hrbuf = malloc(nthreads*nfpert*hlen*sizeof(double));
	TEST_ALLOC(Hrbuf);
	double* const  Hfbuf = malloc(nthreads*nfpert*hlen*sizeof(double));
	TEST_ALLOC(Hfbuf);
	double* const  DDbuf = malloc(nthreads*nfpert*hlen*sizeof(double));
	TEST_ALLOC(DDbuf);

	// allocate buffer for per-filter parameters

	tfarg_t* const tfbuf = malloc(nthreads*nfpert*sizeof(tfarg_t));
	TEST_ALLOC(tfbuf);

	printf("*** Setting up simulation parameters storage\n\n");

	// set up thread arguments

	targ_t targs[nthreads];
	for (size_t i=0; i<nthreads; ++i) {
		targ_t* const targ = &targs[i];

		// thread-independent

		targ->nfpert = nfpert;
		targ->rsize  = rsize;
		targ->fsize  = fsize;
		targ->emmax  = emmax;
		targ->eiff   = eiff;
		targ->tmmax  = tmmax;
		targ->tiff   = tiff;
		targ->tlag   = tlag;

		// thread-dependent

		targ->tnum   = i;
		targ->ebuf   = ebuf  + i*eblen;
		targ->tbuf   = tbuf  + i*tblen;
		targ->tfargs = tfbuf + i*nfpert;

		// thread/filter-dependent

		word_t* const rbufi  = rbuf  + i*nfpert*rlen;
		word_t* const fbufi  = fbuf  + i*nfpert*flen;
		double* const Hrbufi = Hrbuf + i*nfpert*hlen;
		double* const Hfbufi = Hfbuf + i*nfpert*hlen;
		double* const DDbufi = DDbuf + i*nfpert*hlen;
		for (size_t j=0; j<nfpert; ++j) {
			tfarg_t* const tfarg = &targ->tfargs[j];
			rt_randomise(rsize,tfarg->rtab = rbufi+j*rlen,rlam,&rrng);
			rt_randomise(fsize,tfarg->ftab = fbufi+j*flen,flam,&frng);
			tfarg->Hr = Hrbufi+j*hlen;
			tfarg->Hf = Hfbufi+j*hlen;
			tfarg->DD = DDbufi+j*hlen;
		}
	}

	// create threads

	printf("*** Creating %zu threads with %zu simulations per thread\n\n",nthreads,nfpert);

	pthread_t threads[nthreads]; // NOTE: joinable by default (otherwise use pthread_attr_setdetachstate())
	for (size_t i=0; i<nthreads; ++i) {
		const int tres = pthread_create(&threads[i],NULL,compfun,(void*)&targs[i]);
		PASSERT(tres == 0,"unable to create thread %zu",i+1)
	}

	// wait for computational threads to complete

	for (size_t i=0; i<nthreads; ++i) {
		const int tres = pthread_join(threads[i],NULL);
		PASSERT(tres == 0,"unable to join thread %zu",i+1);
	}

	// write out results

	const size_t ofnlen = strlen(odir)+20;
	char ofname[ofnlen];
	snprintf(ofname,ofnlen,"%s/caddr_%zu.dat",odir,jnum);
	printf("\n*** Writing results to \"%s\"... ",ofname);
	fflush(stdout);
	FILE* const dfs = fopen(ofname,"w");
	PASSERT(dfs != NULL,"Failed to open output file \"%s\"\n",ofname);
	fprintf(dfs,"# rule    size    = %2d (lambda  = %8.6f)\n"
	            "# filter  size    = %2d (lambda  = %8.6f)\n"
	            "# entropy seqlen  = %2d (advance = %d)\n"
	            "# dynind  seqlen  = %2d (advance = %d, lag = %d)\n"
	            "# sample  size    = %zu\n\n"
	            ,rsize,rlam,fsize,flam,emmax,eiff,tmmax,tiff,tlag,nthreads*nfpert);
	for (size_t i=0; i<nthreads; ++i) {
		const targ_t* const targ = &targs[i];
		for (size_t j=0; j<targ->nfpert; ++j) {
			const tfarg_t* const tfarg = &targ->tfargs[j];
			fprintf(dfs,"# rule id = ");
			rt_fprint_id(rsize,tfarg->rtab,dfs);
			fprintf(dfs,", filter id = ");
			rt_fprint_id(fsize,tfarg->ftab,dfs);
			fputc('\n',dfs);
			for (int m=0; m<(int)hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,tfarg->Hr[m],tfarg->Hf[m],tfarg->DD[m]);
			fputs("\n",dfs);
		}
	}
	if (fclose(dfs) == -1) PEEXIT("Failed to close output file \"%s\"\n",ofname);
	puts("done\n");

	// free buffers

	printf("*** Freeing memory\n");

	free(tfbuf);
	free(DDbuf);
	free(Hfbuf);
	free(Hrbuf);
	free(tbuf);
	free(ebuf);
	free(fbuf);
	free(rbuf);

	return EXIT_SUCCESS;
}

void* compfun(void* arg)
{
	const double wts = get_wall_time();
	const double cts = get_thread_cpu_time ();

	const pthread_t tpid     = pthread_self();
	const targ_t* const targ = (targ_t*)arg;
	const size_t tnum        = targ->tnum;
	const size_t nfpert      = targ->nfpert;

	printf("thread %2zu (%zu) : %zu filters : STARTED\n",tnum+1,tpid,nfpert);
	fflush(stdout);

	const int rsize = targ->rsize;
	const int fsize = targ->fsize;
	const int emmax = targ->emmax;
	const int eiff  = targ->eiff;
	const int tmmax = targ->tmmax;
	const int tiff  = targ->tiff;
	const int tlag  = targ->tlag;

	uint64_t* const ebuf = targ->ebuf;
	uint64_t* const tbuf = targ->tbuf;

	const int hlen   = (emmax > tmmax ? emmax : tmmax)+1;
	const int rfsize = rsize > fsize ? rsize : fsize;

	for (size_t j=0; j<nfpert; ++j) {

		const tfarg_t* const tfarg = &targ->tfargs[j];

		const word_t* const rtab = tfarg->rtab;
		const word_t* const ftab = tfarg->ftab;
		double*       const Hr   = tfarg->Hr;
		double*       const Hf   = tfarg->Hf;
		double*       const DD   = tfarg->DD;

		for (int m=0; m<hlen; ++m) Hr[m] = NAN;
		for (int m=0; m<hlen; ++m) Hf[m] = NAN;
		for (int m=0; m<hlen; ++m) DD[m] = NAN;
		for (int m=rsize;  m<=emmax; ++m) Hr[m] = rt_entro(rsize,rtab,m,eiff,ebuf)/(double)m;
		for (int m=fsize;  m<=emmax; ++m) Hf[m] = rt_entro(fsize,ftab,m,eiff,ebuf)/(double)m;
		for (int m=rfsize; m<=tmmax; ++m) DD[m] = rt_dd   (rsize,rtab,fsize,ftab,m,tiff,tlag,ebuf,tbuf)/(double)m;

		flockfile(stdout); // prevent another thread butting in!
		printf("\tthread %2zu : filter %2zu of %2zu : rule id = ",tnum+1,j+1,nfpert);
		rt_print_id(rsize,rtab);
		printf(", filter id = ");
		rt_print_id(fsize,ftab);
		printf(" : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f\n",Hr[emmax],Hf[emmax],DD[tmmax]);
		fflush(stdout);
		funlockfile(stdout);
	}

	const double cte = get_thread_cpu_time() - cts;
	const double wte = get_wall_time() - wts;

	printf("thread %2zu : FINISHED (cpu time = %.4f, wall time = %.4f)\n",tnum+1,cte,wte);
	fflush(stdout);

	pthread_exit(NULL);
}
