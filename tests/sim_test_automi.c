#include "word.h"
#include "clap.h"
#include "utils.h"

void mw_automi_ref(const size_t n, const word_t* const w, double* const ami, double* const entro)
{
	const size_t m = n*WBITS;
	const size_t q = m/2; // fine, because WBITS even!
	const double fac = 1.0/(double)m;
	const double p0 = fac*mw_nsetbits(n,w);
	*entro = -xlog2x(p0)-xlog2x(1.0-p0);
	double p[4];
	for (size_t k=1;k<q;++k) {
		int bin[4] = {0};
		for (size_t j=0;j<m;++j) {
			const size_t i   = j+k < m ? j+k : j+k-m; // wrap!
			const size_t jj  = j/WBITS;
			const int    jjj = j%WBITS;
			const size_t ii  = i/WBITS;
			const int    iii = i%WBITS;
			++bin[MIIDX(w[ii],iii,w[jj],jjj)];
		}
		for (int a=0;a<4;++a) p[a] = fac*(double)bin[a];
		ami[k-1] = xlog2x(p[0]) + xlog2x(p[1]) + xlog2x(p[2]) + xlog2x(p[3])
		         - xlog2x(p[0]+p[1]) - xlog2x(p[2]+p[3])
		         - xlog2x(p[0]+p[2]) - xlog2x(p[1]+p[3]);
	}
}

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(n,       size_t,  3,            "number of words");
	CLAP_CARG(wbias,   double,  0.5,          "random word bias");
	CLAP_CARG(seed,    ulong,   0,            "random seed (or 0 for unpredictable)");
	CLAP_CARG(disp,    int,     1,            "display result?");
	CLAP_CARG(logs,    int,     0,            "log scale?");
	puts("---------------------------------------------------------------------------------------\n");

	mt_t rng;
	mt_seed(&rng,seed);
/*
	const word_t wi = wd_random(&rng);
	wd_print(wi); putchar('\n');
	const word_t wj = wd_random(&rng);
	wd_print(wj); putchar('\n');

	const int i = 2, j = 5;
	putchar('\n');
	wd_print(wi>>i); putchar('\n');
	wd_print((wj>>j)<<1); putchar('\n');
	putchar('\n');
	wd_print((wi>>i) & ((word_t)1)); putchar('\n');
	wd_print(((wj>>j)<<1) & ((word_t)2)); putchar('\n');

	size_t idx = MIIDX(wi,i,wj,j);
	printf("\nidx = %zu\n\n",idx);
*/
	const size_t m = WBITS*n;
	const size_t q = m/2; // fine, because WBITS even!

	word_t* const w = mw_alloc(n);
	mw_randomiseb(n,w,wbias,&rng);
	mw_print(n,w);
	putchar('\n');
	putchar('\n');

	double* const amir = calloc(q-1,sizeof(double));
	double* const amio = calloc(q-1,sizeof(double));

	double entr, ento;

	double ts,te;

	ts = timer();
	mw_automi_ref(n,w,amir,&entr);
	te = timer();
	printf("ami ref time = %8.6f\n",te-ts);

	ts = timer();
	mw_automi(n,w,amio,&ento);
	te = timer();
	printf("ami opt time = %8.6f\n",te-ts);

	printf("\nami max abs diff = %.4e (%.4e)\n\n", maxabdiff(q-1,amir,amio),fabs(entr-ento));

for (size_t k=0;k<q-1;++k) printf("%4zu  % 8.6f  % 8.6f\n",k+1,amio[k],ento-amio[k]);

	free(amio);
	free(amir);
	free(w);

	return EXIT_SUCCESS;
}
