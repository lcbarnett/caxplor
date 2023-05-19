#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int sim_viz   (int argc, char* argv[]);
int sim_ana   (int argc, char* argv[]);
int sim_bmark (int argc, char* argv[]);
int sim_test  (int argc, char* argv[]);

int main(int argc, char* argv[])
{
	if (!(argc>1)) {
		fprintf(stderr,"%s: must be at least one program argument\n",argv[0]);
		return EXIT_FAILURE;
	}

	if (strcmp(argv[1],"viz"  )  == 0) return sim_viz   (argc-2,argv+2);
	if (strcmp(argv[1],"ana"  )  == 0) return sim_ana   (argc-2,argv+2);
	if (strcmp(argv[1],"bmark")  == 0) return sim_bmark (argc-2,argv+2);
	if (strcmp(argv[1],"test" )  == 0) return sim_test  (argc-2,argv+2);

	fprintf(stderr,"%s: unknown simulation '%s'\n",argv[0],argv[1]);
	return EXIT_FAILURE;
}
