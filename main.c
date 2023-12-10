#include <stdio.h>
#include <string.h>
#include <stdlib.h>

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

typedef int (*sim_t)(int argc, char* argv[], int info);

int main(int argc, char* argv[])
{
	const double wts = get_wall_time();
	const double cts = get_proc_cpu_time ();

    struct tm loctime = *localtime(&(time_t){time(NULL)});

	if (argc == 1) { // no command line arguments: display compilation options and available simulations
		puts("\ncaxplor compile options:");
#ifdef HAVE_PTHREADS
		puts("\t+WITH_PTHREADS");
#else
		puts("\t-WITH_PTHREADS");
#endif
#ifdef HAVE_X11
		puts("\t+WITH_X11");
#else
		puts("\t-WITH_X11");
#endif
#ifdef HAVE_GD
		puts("\t+WITH_GD");
#else
		puts("\t-WITH_GD");
#endif
		puts("\ncaxplor available simulations:\n\tana\n\tbmark\n\ttest");
#ifdef HAVE_X11
		puts("\txplor");
#endif
#ifdef HAVE_PTHREADS
		puts("\tddf");
		puts("\tddr");
#endif
		putchar('\n');
		return EXIT_SUCCESS;
	}

	// so argc > 1: check simulation and run it

	int res;
	sim_t sim;
	printf("\ncaxplor %s: %s",argv[1],asctime(&loctime));
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
		fprintf(stderr,"%s: unknown simulation \"%s\"\n",argv[0],argv[1]);
		return EXIT_FAILURE;
	}
	if (argc > 2 && strcmp(argv[2],"-i") == 0) { // "dry run": display switches (and optionally other info) and exit
		sim(argc-3,argv+3,1);
		return EXIT_SUCCESS;
	}
	res = sim(argc-2,argv+2,0); // run the sim

	const double cte = get_proc_cpu_time() - cts;
	const double wte = get_wall_time() - wts;

	printf("\nCPU  time (secs) = %.4f\n",cte);
	printf("Wall time (secs) = %.4f\n\n",wte);

	return res;
}
