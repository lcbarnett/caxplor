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
