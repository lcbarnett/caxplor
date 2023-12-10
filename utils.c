#include <string.h>
#include <stdio.h>
#include <time.h>

void report_compilation_options()
{
	puts("\ncaxplor compile options:");
#ifdef HAVE_PTHREADS
	puts("\t+WITH_PTHREADS");
#else
	puts("\t-WITH_PTHREADS");
#endif
#ifdef HAVE_X11
	puts("\t+WITH_X11");
#else
	puts("\t-WITH_X11");
#endif
#ifdef HAVE_GD
	puts("\t+WITH_GD");
#else
	puts("\t-WITH_GD");
#endif
	puts("\ncaxplor available simulations:\n\tana\n\tbmark\n\ttest");
#ifdef HAVE_X11
	puts("\txplor");
#endif
#ifdef HAVE_PTHREADS
	puts("\tddf");
	puts("\tddr");
#endif
	putchar('\n');
}

#ifdef __linux__
	#include <sys/sysinfo.h>
#endif

#ifdef _WIN32
	#include <Windows.h>
#endif

#ifdef __APPLE__
	#include "TargetConditionals.h"
    #if TARGET_IPHONE_SIMULATOR || TARGET_OS_MACCATALYST || TARGET_OS_IPHONE
        // nada
    #elif TARGET_OS_MAC
        #define _OSX
    #else
		#error "Unknown Apple platform"
    #endif
#endif

#include "utils.h"

#define GPDEFTITLE "CA Xplorer"
#define GPCMD "gnuplot -p"
#define SMAXLEN 200

ulong get_free_ram()
{
#if defined __linux__
	struct sysinfo info;
	PASSERT(sysinfo(&info) == 0,"sysinfo call failed");
	return info.freeram;
#elif defined _OSX || defined _WIN32
	return 0; // TODO
#else
	#error Unhandled OS
#endif
}

double get_wall_time()
{
#if defined __linux__ || defined _OSX
	struct timespec ts;
	PASSERT(clock_gettime(CLOCK_REALTIME,&ts) == 0,"'clock_gettime' failed");
	return (double)ts.tv_sec + (double)ts.tv_nsec/1000000000.0;
#elif defined _WIN32
	LARGE_INTEGER time,freq;
	if (!QueryPerformanceFrequency(&freq)) PEEXIT("'QueryPerformanceFrequency' failed");
	if (!QueryPerformanceCounter(&time))   PEEXIT("'QueryPerformanceCounter' failed");
	return (double)time.QuadPart/(double)freq.QuadPart;
#else
	#error Unhandled OS
#endif
}

double get_proc_cpu_time()
{
#if defined __linux__ || defined _OSX
	struct timespec ts;
	PASSERT(clock_gettime(CLOCK_PROCESS_CPUTIME_ID,&ts) == 0,"'clock_gettime' failed");
	return (double)ts.tv_sec + (double)ts.tv_nsec/1000000000.0;
#elif defined _WIN32
	return NAN; // TODO
#else
	#error Unhandled OS
#endif
}

double get_thread_cpu_time()
{
#ifdef __linux__
	struct timespec ts;
	PASSERT(clock_gettime(CLOCK_THREAD_CPUTIME_ID,&ts) == 0,"'clock_gettime' failed");
	return (double)ts.tv_sec + (double)ts.tv_nsec/1000000000.0;
#elif defined _WIN32
	return NAN; // TODO
#else
	#error Unhandled OS
#endif
}

double timer() // deprecate!
{
	return (double)clock()/(double)CLOCKS_PER_SEC;
}

void hist(const size_t n, const double* const x, const size_t m, ulong* const  bin)
{
	const double a = min(n,x);
	const double b = max(n,x);
	const double f = (double)m/(b-a);
	for (size_t i=0;i<m;++i) bin[i] = 0;
	for (size_t k=0;k<n;++k) {
		const ulong j =(ulong)(f*(x[k]-a));
		++bin[j<m?j:m-1];
	}
}

double entro2(const size_t n, const double* const x)
{
	double y = 0.0;
	for (const double* p=x; p<x+n; ++p) y += xlog2x(*p);
	return -y;
}

double* dft_cstab_alloc(const size_t n)
{
	double* const costab = calloc(2*n*n,sizeof(double));
	TEST_ALLOC(costab);
	double* const sintab = costab+n*n; // sin table is offset by n^2
	const double fac = (double)(2.0*M_PI)/(double)n;
	// build sin and cos tables
	for (size_t i=0; i<n; ++i) {
		double* const ctni = costab+n*i;
		double* const stni = sintab+n*i;
		for (size_t j=i; j<n; ++j) sincos(fac*(double)(i*j),stni+j,ctni+j);
	}
	// symmetrise cos table across diagonal
	for (size_t i=0; i<n; ++i) {
		double* const cti  = costab+i;
		double* const ctni = costab+n*i;
		for (size_t j=i+1; j<n; ++j) cti[n*j] = ctni[j];
	}
	// symmetrise sin table across diagonal
	for (size_t i=0; i<n; ++i) {
		double* const sti  = sintab+i;
		double* const stni = sintab+n*i;
		for (size_t j=i+1; j<n; ++j) sti[n*j] = stni[j];
	}
	return costab;
}

void ac2dps(const size_t m, double* const dps, const double* const ac, const double* const costab)
{
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t k=0; k<q; ++k) {
		const size_t mk = m*k;
		double dpsk = 0.0;
		for (size_t j=0; j<q; ++j) dpsk += ac[j  ]*costab[mk+j];
		for (size_t j=q; j<m; ++j) dpsk += ac[m-j]*costab[mk+j];
		dps[k] = dpsk;
	}
}

// statistics

double mean(const size_t n, double* const x, double* const var, const int unbiased)
{
	ASSERT(n > 0,"Empty array!");
	if (n == 1) {
		if (var != NULL) *var = NAN; // undefined
		return *x;
	}
	double sum = 0.0;
	for (size_t i=0;i<n;++i) sum += x[i];
	const double m = sum/(double)n;
	if (var != NULL) {
		double sum = 0.0;
		for (size_t i=0;i<n;++i) {
			double y = x[i]-m;
			sum += y*y;
		}
		*var = sum/(double)(unbiased?n-1:n);
	}
	return m;
}

double median(const size_t n, double* const x, double* const mad, const int unsorted)
{
	// NOTE: changes data in x (unless already sorted and no MAD requested)
	ASSERT(n > 0,"Empty array!");
	if (n == 1) {
		if (mad != NULL) *mad = NAN; // undefined
		return *x;
	}
	if (unsorted) qsort_double(n,x);
	const size_t h = n/2;
	const double med = (2*h == n ? (x[h-1]+x[h])/2.0 : x[h]);
	if (mad != NULL) {
		for (size_t i=0;i<n;++i) x[i] = fabs(x[i]-med);
		*mad = median(n,x,NULL,1);
	}
	return med;
}

void quantiles(const size_t n, double* const x, const size_t Q, double* const q, const int unsorted)
{
	// NOTE: changes data in x (unless already sorted)
	//
	// q must be a double array of size Q+1
	// q[0] is the minimum, q[Q] the maximum of x
	ASSERT(n > 0,"Empty array!");
	if (n <  Q) {for (size_t k=0;k<=Q;++k) q[k] = NAN;  return;} // undefined
	if (unsorted) qsort_double(n,x);
	if (n == Q) {for (size_t k=0;k<=Q;++k) q[k] = x[k]; return;} // trivial
	const double z = (double)n/(double)Q;
	q[0] = x[0];
	for (size_t k=1;k<Q;++k) {
		const double u = (double)k*z-0.5;
		const size_t i = (size_t)floor(u);
		const double v = u-(double)i;
		q[k] = (1.0-v)*x[i]+v*x[i+1];
	}
	q[Q] = x[n-1];
}

size_t FDrule(const size_t n, double* const x, const int unsorted) // Freedman-Diaconis rule for #(histogram bins)
{
	// NOTE: changes data in x (unless already sorted)
	double q[5];
	quantiles(n,x,4,q,unsorted); // quartiles
	return (size_t)ceil((cbrt((double)n)/2.0)*((q[4]-q[0])/(q[3]-q[1])));
}

/*********************************************************************/
/*                      Gnuplot stuff                                */
/*********************************************************************/

FILE* gp_fopen(const char* const gpname, const char* const gpdir, const char* const gpterm, const char* const gptitle, const int xsize, const int ysize)
{
	char gpcmdfile[SMAXLEN];
	snprintf(gpcmdfile,SMAXLEN,"%s/%s.gp",gpdir == NULL ? "/tmp" : gpdir,gpname);
	FILE* const gpc = fopen(gpcmdfile,"w");
	if (gpc == NULL) PEEXIT("failed to open Gnuplot command file \"%s\"\n",gpcmdfile);
	gp_setterm(gpc,gpterm,gptitle,xsize,ysize);
	return gpc;
}

FILE* gp_dopen(const char* const gpname, const char* const gpdir)
{
	char gpdatfile[SMAXLEN];
	snprintf(gpdatfile,SMAXLEN,"%s/%s.dat",gpdir == NULL ? "/tmp" : gpdir,gpname);
	FILE* const gpd = fopen(gpdatfile,"w");
	if (gpd == NULL) PEEXIT("failed to open Gnuplot data file \"%s\"\n",gpdatfile);
	return gpd;
}

void gp_fplot(const char* const gpname, const char* const gpdir)
{
	char gprun[SMAXLEN];
	snprintf(gprun,SMAXLEN,"cd %s && " GPCMD " %s.gp",gpdir == NULL ? "/tmp" : gpdir,gpname);
	if (system(gprun) == -1) PEEXIT("Gnuplot plot command \"%s\" failed\n",gprun);
}

FILE* gp_popen(const char* const gpterm, const char* const gptitle, const int xsize, const int ysize)
{
	FILE* const gpp = popen(GPCMD,"w");
	if (gpp == NULL) PEEXIT("failed to open pipe to Gnuplot\n");
	if (setvbuf(gpp,NULL,_IOLBF,0) != 0) PEEXIT("failed to line-buffer pipe to Gnuplot\n");
	gp_setterm(gpp,gpterm,gptitle,xsize,ysize);
	return gpp;
}

void gp_setterm(FILE* const gp, const char* const gpterm, const char* const gptitle, const int xsize, const int ysize)
{
	// if xsize != 0 and ysize == 0, xsize is a percentage scaling factor
	const int scale = (xsize!=0)&(ysize==0);
	const int xxsize = (xsize == 0 ? 640 : scale ? (int)(6.4*(double)xsize) : xsize);
	const int yysize = (ysize == 0 ? scale ? (int)(4.8*(double)xsize) : 480 : ysize);
	if (gpterm == NULL || strcmp(gpterm,"wxt") == 0) fprintf(gp,"set term \"wxt\" size %d,%d nobackground enhanced title \"%s\" \n",xxsize,yysize,gptitle==NULL?GPDEFTITLE:gptitle);
	else if (strcmp(gpterm,"qt" ) == 0)              fprintf(gp,"set term \"qt\" size %d,%d title \"%s\" enhanced\n",xxsize,yysize,gptitle==NULL?GPDEFTITLE:gptitle);
	else if (strcmp(gpterm,"x11") == 0)              fprintf(gp,"set term \"x11\" title \"%s\" enhanced size %d,%d\n",gptitle==NULL?GPDEFTITLE:gptitle,xxsize,yysize);
	else EEXIT("Unknown Gnuplot terminal \"%s\"",gptitle);
}

void gp_pclose(FILE* const gpp)
{
	if (pclose(gpp) == -1) PEEXIT("failed to close pipe to Gnuplot\n");
}

void gp_binary_write(FILE* const gpp, const size_t n, const double* const x, const int inplace)
{
	// WARNING: if 'inplace' is set, x is unusable after calling (but should still be freed if allocated on heap)!
	if (inplace) {
		const float* const xf = double2float_inplace(n,x);
		fwrite(xf,sizeof(float),n,gpp);
	}
	else {
		float* const xf = double2float_alloc(n,x);
		fwrite(xf,sizeof(float),n,gpp);
		free(xf);
	}
}

const char* gp_palette[] = { // see also https://github.com/Gnuplotting/gnuplot-palettes
	"0 '#000090', 1 '#000fff', 2 '#0090ff', 3 '#0fffee',  4 '#90ff70',   5 '#ffee00',    6 '#ff7000', 7 '#ee0000', 8 '#7f0000'", // Matlab-like :-)
	"0 '#FFF7EC', 1 '#FEE8C8', 2 '#FDD49E', 3 '#FDBB84',  4 '#FC8D59',   5 '#EF6548',    6 '#D7301F', 7 '#990000'", // OrRd
	"0 'black',   3 'blue',    6 'green',   9 'yellow',  12 'orange',   15 'red',      100 'dark-red'",
	"0 'black',   1 'blue',    3 'green',   6 'yellow',  10 'orange',   15 'red',      100 'dark-red'",
	"0 '#00008f', 2 '#0000ff', 4 '#00ffff', 8 '#ffff00', 24 '#ff0000', 128 '#800000'"
};

#undef GPDEFTITLE
#undef GPCMD
#undef SMAXLEN
