#include "word.h"
#include "clap.h"
#include "utils.h"

int sim_test(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(ncols,   int,     10,           "number of columns");
	CLAP_CARG(nrows,   int,     20,           "number of rows");
	puts("---------------------------------------------------------------------------------------\n");

	FILE* data = gp_dopen("test",NULL);
	for (int i=0;i<nrows;++i) {
		for (int j=0;j<ncols;++j) {
			const float p = sinf((float)i)*cosf((float)j);
			fwrite(&p,sizeof(float),1,data);
		}
	}
	fclose(data);

	FILE* gp = gp_fopen("test",NULL,NULL);
	fprintf(gp,"data = \"test.dat\"\n");
	fprintf(gp,"set xr [-0.5:%g]\n",(double)ncols-0.5);
	fprintf(gp,"set yr [-0.5:%g]\n",(double)nrows-0.5);
	fprintf(gp,"plot data binary array=(%d,%d) with image\n",ncols,nrows);
	fclose(gp);

	gp_fplot("test",NULL,NULL);

	return EXIT_SUCCESS;
}
