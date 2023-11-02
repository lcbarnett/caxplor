#ifndef CAUTILS
#define CAUTILS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

/*********************************************************************/
/*                      Useful macros                                */
/*********************************************************************/

// bomb out gracefully with diagnostics

#define ERRPT fprintf(stderr,"ERROR in '%s' [%s:%d]: ",__FUNCTION__,__FILE__,__LINE__)

#define EEXIT(...)  {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}
#define PEEXIT(...) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}

#define ASSERT(cond,...)  {if (!(cond)) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}}
#define PASSERT(cond,...) {if (!(cond)) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}}

/*********************************************************************/

// entropy stuff

static inline double xlog2x(const double x)
{
	return (x > DBL_MIN ? x*log2(x) : (x > -DBL_MIN ? 0.0 : NAN));
}

double entro2(const size_t n, const double* const x);

// cosine and sin tables for DFT

double* dft_cstab_alloc (const size_t n); // remember to free return!

void ac2dps(const size_t n, double* const dps, const double* const ac, const double* const costab);

// sort scalar arrays

#define QSORT_COMP(type) \
static inline int qsort_##type##_comp(const void* const x1, const void* const x2) \
{ \
    return (*(const type* const)x2 < *(const type* const)x1 ? 1 : *(const type* const)x2 > *(const type* const)x1 ? -1 : 0); \
}

#define QSORT(type) \
static inline void qsort_##type(const size_t n, type* const x) \
{ \
	qsort(x,n,sizeof(type),qsort_##type##_comp); \
}

#define QSORT_DEFINE(type) QSORT_COMP(type) QSORT(type)

// instantiate the ones you want

QSORT_DEFINE(double)
QSORT_DEFINE(int)

void hist(const size_t n, const double* const x, const size_t m, ulong* const  bin);

static inline float* double2float_alloc(const size_t n, const double* const x)
{
	float* const xf = malloc(n*sizeof(float)); // !!! remember to free !!!
	float* pf = xf;
	for (const double* p=x;p<x+n;++p,++pf) *pf = (float)*p;
	return xf;
}

static inline float* double2float_inplace(const size_t n, const double* const x) // in-place version
{
	// WARNING: This code is mad, bad and dangerous: note that it leaves x unusable
	float* const xf = (float*)x; // alias as array of float
	float* pf = xf;
	for (const double* p=x;p<x+n;++p,++pf) *pf = (float)*p; // "okay" because sizeof(float) <= sizeof(double)
	return xf;
}

static inline double max(const size_t n, const double* const x)
{
	double d = -INFINITY;
	for (size_t i=0;i<n;++i) if (x[i] > d) d = x[i];
	return d;
}

static inline double min(const size_t n, const double* const x)
{
	double d = INFINITY;
	for (size_t i=0;i<n;++i) if (x[i] < d) d = x[i];
	return d;
}

static inline double maxabs(const size_t n, const double* const x)
{
	double d = 0.0;
	for (size_t i=0;i<n;++i) {
		const double di = fabs(x[i]);
		if (di > d) d = di;
	}
	return d;
}

static inline double maxabdiff(const size_t n, const double* const x, const double* const y)
{
	double d = 0.0;
	for (size_t i=0;i<n;++i) {
		const double di = fabs(x[i]-y[i]);
		if (di > d) d = di;
	}
	return d;
}

static inline float maxf(const size_t n, const float* const x)
{
	float d = -INFINITY;
	for (size_t i=0;i<n;++i) if (x[i] > d) d = x[i];
	return d;
}

static inline float minf(const size_t n, const float* const x)
{
	float d = INFINITY;
	for (size_t i=0;i<n;++i) if (x[i] < d) d = x[i];
	return d;
}

static inline float maxabsf(const size_t n, const float* const x)
{
	float d = 0.0;
	for (size_t i=0;i<n;++i) {
		const float di = fabsf(x[i]);
		if (di > d) d = di;
	}
	return d;
}

static inline float maxabdifff(const size_t n, const float* const x, const float* const y)
{
	float d = 0.0;
	for (size_t i=0;i<n;++i) {
		const float di = fabsf(x[i]-y[i]);
		if (di > d) d = di;
	}
	return d;
}

static inline void scale(const size_t n, double* const x, const double fac)
{
	for (double* p=x;p<x+n;++p) *p *= fac;
}

static inline void scalef(const size_t n, float* const x, const float fac)
{
	for (float*  p=x;p<x+n;++p) *p *= fac;
}

static inline void sqmag(const size_t n, double* const a, const double* const x, const double* const y)
{
	for (size_t i=0;i<n;++i) a[i] = x[i]*x[i]+y[i]*y[i];
}

static inline void sqmagf(const size_t n, float* const a, const float* const x, const float* const y)
{
	for (size_t i=0;i<n;++i) a[i] = x[i]*x[i]+y[i]*y[i];
}

double get_wall_time();
double get_cpu_time();
double timer();

// statistics

double mean      (const size_t n, double* const x, double* const var, const int unbiased);
double median    (const size_t n, double* const x, double* const mad, const int unsorted);
void   quantiles (const size_t n, double* const x, const size_t Q, double* const q, const int unsorted);
size_t FDrule    (const size_t n, double* const x, const int unsorted);

/*********************************************************************/
/*                      Gnuplot stuff                                */
/*********************************************************************/

FILE* gp_dopen(const char* const gpname, const char* const gpdir);
FILE* gp_fopen(const char* const gpname, const char* const gpdir, const char* const gpterm, const char* const gptitle, const int xsize, const int ysize);
FILE* gp_popen(const char* const gpterm, const char* const gptitle, const int xsize, const int ysize);
void  gp_pclose(FILE* const gpp);
void  gp_fplot(const char* const gpname, const char* const gpdir);
void  gp_setterm(FILE* const gp, const char* const gpterm, const char* const gptitle, const int xsize, const int ysize);


// WARNING: if 'inplace' is set, x is unusable after calling (but should still be freed if allocated on heap)!
void  gp_binary_write(FILE* const gpp, const size_t n, const double* const x, const int inplace);

extern const char* gp_palette[];

#endif // CAUTILS
