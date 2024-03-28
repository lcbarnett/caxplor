#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "utils.h"

int sim_ana   (int argc, char* argv[], int info);
int sim_bmark (int argc, char* argv[], int info);
int sim_test  (int argc, char* argv[], int info);
#ifdef HAVE_X11
int sim_xplor (int argc, char* argv[], int info);
#endif
#ifdef HAVE_PTHREADS
int sim_ddf   (int argc, char* argv[], int info);
int sim_ddr   (int argc, char* argv[], int info);
#endif

int main(int argc, char* argv[])
{
	// if no command line arguments display compilation options and available simulations and exit

	if (argc == 1) {
		report_compilation_options();
		return EXIT_SUCCESS;
	}

	// timings

	const double wts = get_wall_time();
	const double cts = get_proc_cpu_time ();

	printf("\ncaxplor %s: %s",argv[1],asctime(localtime(&(time_t){time(NULL)})));

	// find which simulation requested (note argc > 1)

	int (*sim)(int argc, char* argv[], int info); // pointer to simulation function

	if      (strcmp(argv[1],"ana"  )  == 0) sim = sim_ana;
	else if (strcmp(argv[1],"bmark")  == 0) sim = sim_bmark;
	else if (strcmp(argv[1],"test" )  == 0) sim = sim_test;
#ifdef HAVE_X11
	else if (strcmp(argv[1],"xplor")  == 0) sim = sim_xplor;
#endif
#ifdef HAVE_PTHREADS
	else if (strcmp(argv[1],"ddf"  )  == 0) sim = sim_ddf;
	else if (strcmp(argv[1],"ddr"  )  == 0) sim = sim_ddr;
#endif
	else {
		fprintf(stderr,"\nUnknown simulation \"%s\"\n",argv[1]);
		return EXIT_FAILURE;
	}

	// simulation found; first switch -i specifies "dry run": display parameters (and optionally other info) and exit

	if (argc > 2 && strcmp(argv[2],"-i") == 0) {
		sim(argc-3,argv+3,1); // run with info  == 1 (true)
		return EXIT_SUCCESS;
	}

	// run the simulation

	const int res = sim(argc-2,argv+2,0);

	// report timings

	const double cte = get_proc_cpu_time() - cts;
	const double wte = get_wall_time() - wts;

	printf("\nCPU  time (secs) = %.4f\n",cte);
	printf("Wall time (secs) = %.4f\n\n",wte);

	return res;
}
