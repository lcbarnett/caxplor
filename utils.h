#ifndef CAUTILS
#define CAUTILS

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

// entropy stuff

static inline double xlog2x(const double x)
{
	return (x > DBL_MIN ? x*log2(x) : (x > -DBL_MIN ? 0.0 : NAN));
}

double entro2(const size_t n, const double* const x);

// cosine and sin tables for DFT - remember to free return!

double* dft_costab_alloc(const size_t n);
double* dft_sintab_alloc(const size_t n);

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
