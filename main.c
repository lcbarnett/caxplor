#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "utils.h"

int sim_ana   (int argc, char* argv[]);
int sim_bmark (int argc, char* argv[]);
int sim_test  (int argc, char* argv[]);
#ifdef HAVE_X11
int sim_xplor (int argc, char* argv[]);
#endif
#ifdef HAVE_PTHREADS
int sim_dd_mt (int argc, char* argv[]);
#endif

int main(int argc, char* argv[])
{
	const double wts = get_wall_time();
	const double cts = get_cpu_time ();

	int res;
	if (argc > 1) {
		if      (strcmp(argv[1],"ana"  )  == 0) res = sim_ana   (argc-2,argv+2);
		else if (strcmp(argv[1],"bmark")  == 0) res = sim_bmark (argc-2,argv+2);
		else if (strcmp(argv[1],"test" )  == 0) res = sim_test  (argc-2,argv+2);
#ifdef HAVE_X11
		else if (strcmp(argv[1],"xplor")  == 0) res = sim_xplor (argc-2,argv+2);
#endif
#ifdef HAVE_PTHREADS
		else if (strcmp(argv[1],"dd"   )  == 0) res = sim_dd_mt (argc-2,argv+2);
#endif
		else {
			fprintf(stderr,"%s: unknown simulation \"%s\"\n",argv[0],argv[1]);
			return EXIT_FAILURE;
		}
	}
	else {
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
		puts("\tdd");
#endif
		putchar('\n');
		return EXIT_SUCCESS;
	}

	const double cte = get_cpu_time () - cts;
	const double wte = get_wall_time() - wts;

	printf("\nCPU  time (secs) = %.4f\n",cte);
	printf("Wall time (secs) = %.4f\n\n",wte);

	return res;
}
