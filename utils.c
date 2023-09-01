#include <string.h>

#include "utils.h"

#define GPDEFTERM "wxt size 640,480 nobackground enhanced title 'CA Explorer' persist raise"
#define SMAXLEN 200

double entro2(const size_t n, const double* const x)
{
	double y = 0.0;
	for (const double* p=x; p<x+n; ++p) y += xlog2x(*p);
	return -y;
}

double* dft_costab_alloc(const size_t n)
{
	double* const costab = calloc(n*n,sizeof(double));
	const double fac = (2.0*M_PI)/(double)n;
	for (size_t i=0; i<n; ++i) {
		double* const ctni = costab+n*i;
		for (size_t j=i; j<n; ++j) ctni[j] = cos(fac*(double)(i*j));
	}
	// symmetrise across diagonal
	for (size_t i=0; i<n; ++i) {
		double* const cti  = costab+i;
		double* const ctni = costab+n*i;
		for (size_t j=i+1; j<n; ++j) cti[n*j] = ctni[j];
	}
	return costab;
}

double* dft_sintab_alloc(const size_t n)
{
	double* const sintab = calloc(n*n,sizeof(double));
	const double fac = (2.0*M_PI)/(double)n;
	for (size_t i=0; i<n; ++i) {
		double* const stni = sintab+n*i;
		for (size_t j=i; j<n; ++j) stni[j] = sin(fac*(double)(i*j));
	}
	// symmetrise across diagonal
	for (size_t i=0; i<n; ++i) {
		double* const sti  = sintab+i;
		double* const stni = sintab+n*i;
		for (size_t j=i+1; j<n; ++j) sti[n*j] = stni[j];
	}
	return sintab;
}

void dps(const size_t n, double* const s, const double* const ac, const double* const costab)
{
	//const double ac0 = ac[0];
	for (size_t k=0; k<n; ++k) {
		const size_t nk = n*k;
		double sk = 0.0;
		for (size_t j=0; j<n; ++j) sk += ac[j]*costab[nk+j];
		//s[k] = 2.0*sk-ac0;
		s[k] = sk;
	}
}

/*********************************************************************/
/*                      Gnuplot stuff                                */
/*********************************************************************/

FILE* gp_fopen(const char* const gpname, const char* const gpdir, const char* const gpterm)
{
	char gpcmdfile[SMAXLEN];
	snprintf(gpcmdfile,SMAXLEN,"%s/%s.gp",gpdir == NULL ? "/tmp" : gpdir,gpname);
	FILE* const gpc = fopen(gpcmdfile,"w");
	if (gpc == NULL) PEEXIT("failed to open Gnuplot command file \"%s\"\n",gpcmdfile);
	fprintf(gpc,"set term %s\n",gpterm == NULL ? GPDEFTERM : gpterm);
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

void gp_fplot(const char* const gpname, const char* const gpdir, const char* const gpcmd)
{
	char gprun[SMAXLEN];
	snprintf(gprun,SMAXLEN,"cd %s && %s %s.gp",gpdir == NULL ? "/tmp" : gpdir,gpcmd == NULL ? "gnuplot -p" : gpcmd,gpname);
	if (system(gprun) == -1) PEEXIT("Gnuplot plot command \"%s\" failed\n",gprun);
}

FILE* gp_popen(const char* const gpcmd, const char* const gpterm)
{
	FILE* const gpp = popen(gpcmd == NULL ? "gnuplot" : gpcmd,"w");
	if (gpp == NULL) PEEXIT("failed to open pipe to Gnuplot\n");
	if (setvbuf(gpp,NULL,_IOLBF,0) != 0) PEEXIT("failed to line-buffer pipe to Gnuplot\n");
	fprintf(gpp,"set term %s\n",gpterm == NULL ? GPDEFTERM : gpterm);
	return gpp;
}

void gp_pclose(FILE* const gpp)
{
	if (pclose(gpp) == -1) PEEXIT("failed to close pipe to Gnuplot\n");
}

#undef GPDEFTERM
