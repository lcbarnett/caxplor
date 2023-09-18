#include <string.h>

#include "utils.h"

#define GPDEFTERM "wxt size 640,480 nobackground enhanced title 'Gnuplot: CA Xplorer' persist raise"
#define SMAXLEN 200

double entro2(const size_t n, const double* const x)
{
	double y = 0.0;
	for (const double* p=x; p<x+n; ++p) y += xlog2x(*p);
	return -y;
}

double* dft_cstab_alloc(const size_t n)
{
	double* const costab = calloc(2*n*n,sizeof(double));
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
//fprintf(stderr,"q = %4zu  k = %4zu\m",q,k);
		dps[k] = dpsk;
	}
//for (size_t k=0;k<q;++k) printf("%4zu  % 16.4f  % 16.4f\n",k,ac[k],dps[k]);
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

const char* gp_palette[] = {
	"0 'black', 3 'blue', 6 'green', 9 'yellow', 12 'orange', 15 'red', 100 'dark-red'",
	"0 'black', 1 'blue', 3 'green', 6 'yellow', 10 'orange', 15 'red', 100 'dark-red'",
	"1 '#00008f', 8 '#0000ff', 24 '#00ffff', 40 '#ffff00', 56 '#ff0000', 64 '#800000'"
};

#undef GPDEFTERM
