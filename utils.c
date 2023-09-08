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

dft_float_t* dft_cstab_alloc(const size_t n)
{
	dft_float_t* const costab = calloc(2*n*n,sizeof(dft_float_t));
	dft_float_t* const sintab = costab+n*n; // sin table is offset by n^2
	const dft_float_t fac = (dft_float_t)(2.0*M_PI)/(dft_float_t)n;
	// build sin and cos tables
	for (size_t i=0; i<n; ++i) {
		dft_float_t* const ctni = costab+n*i;
		dft_float_t* const stni = sintab+n*i;
		for (size_t j=i; j<n; ++j) {
#ifdef DFT_SINGLE_PREC_FLOAT
			sincosf(fac*(dft_float_t)(i*j),stni+j,ctni+j);
#else
			sincos(fac*(dft_float_t)(i*j),stni+j,ctni+j);
#endif
		}
	}
	// symmetrise cos table across diagonal
	for (size_t i=0; i<n; ++i) {
		dft_float_t* const cti  = costab+i;
		dft_float_t* const ctni = costab+n*i;
		for (size_t j=i+1; j<n; ++j) cti[n*j] = ctni[j];
	}
	// symmetrise sin table across diagonal
	for (size_t i=0; i<n; ++i) {
		dft_float_t* const sti  = sintab+i;
		dft_float_t* const stni = sintab+n*i;
		for (size_t j=i+1; j<n; ++j) sti[n*j] = stni[j];
	}
	return costab;
}

void ac2dps(const size_t n, dft_float_t* const dps, const dft_float_t* const ac, const dft_float_t* const costab)
{
	//const double ac0 = ac[0];
	for (size_t k=0; k<n; ++k) {
		const size_t nk = n*k;
		dft_float_t sk = 0.0;
		for (size_t j=0; j<n; ++j) sk += ac[j]*costab[nk+j];
		//dps[k] = 2.0*sk-ac0;
		dps[k] = sk;
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
