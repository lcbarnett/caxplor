#include "analyse.h"
#include "utils.h"
#include "ca.h"

void caana_automi
(
	// sizes
	const size_t        n,
	const size_t        I,
	const size_t        q,
	const size_t        Q,
	// CAs
	const word_t* const ca,
	const word_t* const fca,
	// flags
	const int           filtering,
	const int           amice,
	// Gnuplot parameters
	const int           gpipw
)
{
	printf("calculating CA auto-MI ... "); fflush(stdout);

	// calculate aut-MI
	double* const ami = malloc(Q*sizeof(double));
	if (filtering) ca_automi(I,n,fca,ami); else ca_automi(I,n,ca,ami);
	if (amice) {
		for (size_t r=0;r<Q;r+=q) {
			for (size_t k=1;k<q;++k) ami[r+k] -= ami[r];
		}
	}
	for (size_t i=0;i<Q;i+=q) ami[i] = NAN; // suppress I(0)
	const double amimax = max(Q,ami);
	const double amimin = min(Q,ami);

	// calculate medians per distance
	double* const amik = malloc(I*sizeof(double));
	double* const med  = malloc(q*sizeof(double));
	med[0] = NAN;
	for (size_t k=1;k<q;++k) {
		for (size_t r=0;r<I;++r) amik[r] = ami[q*r+k];
		med[k] = median(I,amik,NULL,0);
	}

	// display AMI heat map
	printf("AMI max = %g, min = %g\n",amimax,amimin);
	FILE* const gp1 = gp_popen(NULL,NULL);
	fprintf(gp1,"set size ratio -1\n");
	fprintf(gp1,"unset xtics\n");
	fprintf(gp1,"unset ytics\n");
	fprintf(gp1,"set palette defined (%s)\n",gp_palette[2]);
	fprintf(gp1,"set cbr [0:*]\n");
	fprintf(gp1,"set xr [+0.5:%g]\n",(double)q-0.5);
	fprintf(gp1,"set yr [-0.5:%g]\n",(double)I-0.5);
	fprintf(gp1,"plot '-' binary array=(%zu,%zu) flip=y with image not\n",q,I);
	gp_binary_write(gp1,Q,ami,gpipw); // NOTE: if gpipw set, ami is now unusable!
	if (pclose(gp1) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");

	// display medians per distance box plot
	FILE* const gp2 = gp_popen(NULL,NULL);
	fprintf(gp2,"set xr [0.5:%g]\n",(double)q+0.5);
	fprintf(gp2,"set style fill solid 1\n");
	fprintf(gp2,"plot '-' u 1:2 w boxes not\n");
	for (size_t k=1;k<q;++k) fprintf(gp2,"%4zu  %12.8f f\n",k,med[k]);
	if (pclose(gp2) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");

	free(med);
	free(amik);
	free(ami);
}
