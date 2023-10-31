#include <pthread.h>

#include "clap.h"
#include "rtab.h"

#define ofnlen 1000

typedef struct {
	size_t tnum;
	rtl_t* rule;
	rtl_t* filt;
	int    emmax;
	int    eiff;
	int    tmmax;
	int    tiff;
	int    tlag;
	char   ofname[ofnlen+1];
} targ_t;

void* compfun(void* arg)
{
	const pthread_t tpid = pthread_self();
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
	const char* const ofname = targ->ofname;

	printf("thread %2zu (%zu): rule id = ",tnum,tpid);
	rt_print_id(rsize,rtab);
	printf(", filter id = ");
	rt_print_id(fsize,ftab);
	putchar('\n'); fflush(stdout);

	const int hlen   = (emmax > tmmax ? emmax : tmmax)+1;
	const int rfsize = rsize > fsize ? rsize : fsize;

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

	FILE* const dfs = fopen(ofname,"w");
	if (dfs == NULL) PEEXIT("thread %2zu : failed to open output file \"%s\"\n",tnum,ofname);
	for (int m=0; m<hlen; ++m) fprintf(dfs,"%4d\t%8.6f\t%8.6f\t%8.6f\n",m,Hr[m],Hf[m],DD[m]);
	if (fclose(dfs) == -1) PEEXIT("thread %2zu : failed to close output file\n",tnum);

	printf("thread %2zu : rule entropy ≈ %8.6f, filter entropy ≈ %8.6f, DD ≈ %8.6f : written \"%s\"\n",tnum,Hr[emmax],Hf[emmax],DD[tmmax],ofname);

	pthread_exit(NULL);
}
