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
	size_t   tnum;
	size_t   nfpert;
	int      rsize;
	int      fsize;
	int      emmax;
	int      eiff;
	int      tmmax;
	int      tiff;
	int      tlag;
	tfarg_t* tfargs;
} targ_t;

static void* compfun(void* arg);

int sim_ddr(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name      type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_VARG(rsize,    int,     5,             "CA rule size");
	CLAP_VARG(rlam,     double,  0.6,           "CA rule lambda");
	CLAP_CARG(rseed,    ulong,   0,             "CA rule random seed (or 0 for unpredictable)");
	CLAP_VARG(fsize,    int,     5,             "filter rule size");
	CLAP_VARG(flammin,  double,  0.6,           "filter rule lambda range minimum");
	CLAP_VARG(flammax,  double,  0.9,           "filter rule lambda range maximum");
	CLAP_VARG(flamres,  size_t,  10,            "filter rule lambda range resolution");
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

	const char* const jobidxs = getenv(jobidx);
	ASSERT(jobidxs != NULL,"Failed to find environmental variable \"%s\"",jobidx);
	const size_t jnum = (size_t)atoi(jobidxs);
	ASSERT(jnum > 0 && jnum <= flamres, "Bad job index (%zu)",jnum);
	const double flam = flammin+(double)(jnum-1)*((flammax-flammin)/((double)(flamres-1)));
	printf("*** Job number %zu: filter lambda = %g ***\n\n",jnum,flam);

	// pseudo-random number generators

	mt_t rrng, frng;
	mt_seed(&rrng,rseed);
	mt_seed(&frng,fseed);

	// allocate buffers for random rules and filters

	const   size_t rlen = POW2(rsize);
	word_t* const  rbuf = malloc(nthreads*nfpert*rlen*sizeof(word_t));

	const   size_t flen = POW2(fsize);
	word_t* const  fbuf = malloc(nthreads*nfpert*flen*sizeof(word_t));

	// allocate buffers for entropies and DD

	const   size_t hlen  = (size_t)(emmax > tmmax ? emmax : tmmax)+1;
	double* const  Hrbuf = malloc(nthreads*nfpert*hlen*sizeof(double));
	double* const  Hfbuf = malloc(nthreads*nfpert*hlen*sizeof(double));
	double* const  DDbuf = malloc(nthreads*nfpert*hlen*sizeof(double));

	// allocate buffer for per-filter parameters

	const size_t   tflen = sizeof(tfarg_t);
	tfarg_t* const tfbuf = malloc(nthreads*nfpert*tflen);

	// set up thread-independent parameters

	targ_t targs[nthreads];
	for (size_t i=0; i<nthreads; ++i) {
		targs[i].rsize  = rsize;
		targs[i].fsize  = fsize;
		targs[i].emmax  = emmax;
		targs[i].eiff   = eiff;
		targs[i].tmmax  = tmmax;
		targs[i].tiff   = tiff;
		targs[i].tlag   = tlag;
		targs[i].tfargs = tfbuf + i*nfpert*tflen;
		TEST_ALLOC(targs[i].tfargs);
	}

	// loop through rules/filters, setting up thread-dependent parameters

	for (size_t i=0; i<nthreads; ++i) {
		targ_t* const targ = &targs[i];
		targ->tnum   = i;
		targ->nfpert = nfpert;
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
	printf("\nWriting results to \"%s\"... ",ofname);
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
	puts("done");

	// free buffers

	free(tfbuf);
	free(DDbuf);
	free(Hfbuf);
	free(Hrbuf);
	free(fbuf);
	free(rbuf);

	return EXIT_SUCCESS;
}

void* compfun(void* arg)
{
	const pthread_t tpid = pthread_self();
	const targ_t* const targ = (targ_t*)arg;
	const size_t tnum   = targ->tnum;
	const size_t nfpert = targ->nfpert;

	printf("thread %2zu (%zu) : %zu filters : STARTED\n",tnum+1,tpid,nfpert);
	fflush(stdout);

	const int rsize = targ->rsize;
	const int fsize = targ->fsize;
	const int emmax = targ->emmax;
	const int eiff  = targ->eiff;
	const int tmmax = targ->tmmax;
	const int tiff  = targ->tiff;
	const int tlag  = targ->tlag;

	const int hlen   = (emmax > tmmax ? emmax : tmmax)+1;
	const int rfsize = rsize > fsize ? rsize : fsize;

	for (size_t j=0; j<nfpert; ++j) {

		const tfarg_t* const tfarg = &targ->tfargs[j];

		const word_t* const rtab  = tfarg->rtab;
		const word_t* const ftab  = tfarg->ftab;
		double*       const Hr    = tfarg->Hr;
		double*       const Hf    = tfarg->Hf;
		double*       const DD    = tfarg->DD;

		for (int m=0; m<hlen; ++m) Hr[m] = NAN;
		for (int m=0; m<hlen; ++m) Hf[m] = NAN;
		for (int m=0; m<hlen; ++m) DD[m] = NAN;
		for (int m=rsize;  m<=emmax; ++m) Hr[m] = rt_entro(rsize,rtab,m,eiff)/(double)m;
		for (int m=fsize;  m<=emmax; ++m) Hf[m] = rt_entro(fsize,ftab,m,eiff)/(double)m;
		for (int m=rfsize; m<=tmmax; ++m) DD[m] = rt_trent1(rsize,rtab,fsize,ftab,m,tiff,tlag)/(double)m;

		flockfile(stdout); // prevent another thread butting in!
		printf("\tthread %2zu : filter %2zu of %2zu : rule id = ",tnum+1,j+1,nfpert);
		rt_print_id(rsize,rtab);
		printf(", filter id = ");
		rt_print_id(fsize,ftab);
		printf(" : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f\n",Hr[emmax],Hf[emmax],DD[tmmax]);
		fflush(stdout);
		funlockfile(stdout);
	}

	printf("thread %2zu : FINISHED\n",tnum+1);
	fflush(stdout);

	pthread_exit(NULL);
}
