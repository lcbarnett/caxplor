#ifndef CAUTILS
#define CAUTILS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>

// entropy stuff

static inline double xlog2x(const double x)
{
	return (x > DBL_MIN ? x*log2(x) : (x > -DBL_MIN ? 0.0 : NAN));
}

double entro2(const size_t n, const double* const x);

// cosine and sin tables for DFT

#ifdef DFT_SINGLE_PREC_FLOAT
	typedef float  dft_float_t;
#else
	typedef double dft_float_t;
#endif

dft_float_t* dft_cstab_alloc (const size_t n); // remember to free return!

void ac2dps(const size_t n, dft_float_t* const dps, const dft_float_t* const ac, const dft_float_t* const costab);

// misc stuff

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

static inline double timer()
{
	return (double)clock()/(double)CLOCKS_PER_SEC;
}

/*********************************************************************/
/*                      Gnuplot stuff                                */
/*********************************************************************/

FILE* gp_fopen(const char* const gpname, const char* const gpdir, const char* const gpterm);
FILE* gp_dopen(const char* const gpname, const char* const gpdir);
void  gp_fplot(const char* const gpname, const char* const gpdir, const char* const gpcmd);

FILE* gp_popen(const char* const gpcmd,  const char* const gpterm);
void  gp_pclose(FILE* const gpp);


/*********************************************************************/
/*                      Useful macros                                */
/*********************************************************************/

// bomb out gracefully with diagnostics

#define ERRPT fprintf(stderr,"ERROR in '%s' [%s:%d]: ",__FUNCTION__,__FILE__,__LINE__)

#define EEXIT(...)  {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}
#define PEEXIT(...) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}

#define ASSERT(cond,...)  {if (!(cond)) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); exit(EXIT_FAILURE);}}
#define PASSERT(cond,...) {if (!(cond)) {ERRPT; fprintf(stderr,__VA_ARGS__); fputc('\n',stderr); perror(NULL); exit(EXIT_FAILURE);}}

#endif // CAUTILS
