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
	CLAP_CARG(ncols,   int,     100,          "number of columns");
	CLAP_CARG(nrows,   int,     200,          "number of rows");
	CLAP_CARG(gpipe,   int,     0,            "pipe to Gnuplot?");
	puts("---------------------------------------------------------------------------------------\n");

	// The data

	const size_t npix = (size_t)ncols*(size_t)nrows;
	float* const dat = malloc(npix*sizeof(float));
	float* pdat = dat;
	for (int i=0;i<nrows;++i) {
		const float y = (float)i/(float)(nrows-1);
		for (int j=0;j<ncols;++j) {
			const float x = (float)(j-ncols)/(float)(ncols-1);
			*pdat++ = x*cosf(x)*sinf(1.0f/(0.1f+y/10.0f));
		}
	}

	// Note: for Gnuplot binary array with image, default graph 0,0 is bottom left; flip=y for top left

	if (gpipe) {
		printf("Plotting to Gnuplot pipe\n");

		// Pipe data to Gnuplot (binary) and plot as image

		FILE* gp = gp_popen(NULL,NULL);
		fprintf(gp,"set size ratio -1\n");
		fprintf(gp,"set xr [-0.5:%g]\n",(double)ncols-0.5);
		fprintf(gp,"set yr [-0.5:%g]\n",(double)nrows-0.5);
		fprintf(gp,"plot '-' binary array=(%d,%d) flip=y with image not\n",ncols,nrows);
		fwrite(dat,sizeof(float),npix,gp);
		if (pclose(gp) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");
	}
	else {
		printf("Plotting using Gnuplot command file\n");

		// Write data to file (binary)

		FILE* data = gp_dopen("test",NULL);
		fwrite(dat,sizeof(float),npix,data);
		if (fclose(data) == EOF) PEEXIT("failed to close Gnuplot data file\n");

		// Plot data as image

		FILE* gp = gp_fopen("test",NULL,NULL);
		fprintf(gp,"data = \"test.dat\"\n");
		fprintf(gp,"set size ratio -1\n");
		fprintf(gp,"set xr [-0.5:%g]\n",(double)ncols-0.5);
		fprintf(gp,"set yr [-0.5:%g]\n",(double)nrows-0.5);
		fprintf(gp,"plot data binary array=(%d,%d) flip=y with image not\n",ncols,nrows);
		if (fclose(gp) == EOF) PEEXIT("failed to close Gnuplot command file\n");
		gp_fplot("test",NULL,NULL);
	}

	free(dat);

	return EXIT_SUCCESS;
}
