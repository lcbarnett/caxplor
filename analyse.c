#include "analyse.h"
#include "utils.h"
#include "ca.h"

void caana_period
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const rtl_t*  const rule,
	const size_t        prff,
	const size_t        pmax
)
{
	printf("calculating CA period... "); fflush(stdout);
	word_t* const wca = mw_copy_alloc(I*n,ca); // copy to working CA
	mw_run(prff,n,wca,rule->size,rule->tab);
	int prot;
	const size_t period = ca_period(pmax,n,wca,rule->size,rule->tab,&prot);
	if (period == pmax) printf("> %zu iterations\n",pmax); else printf("%zu iterations, twist = %d\n",period,prot);
	free(wca);
}

void caana_dps
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const word_t* const fca,
	const int           filtering,
	const double* const costab,
	const int           gpipw
)
{
	printf("calculating CA spectrum ... "); fflush(stdout);
	const size_t m = n*WBITS; // bits in a CA row
	const size_t q = m/2+1;   // half+1 bits in a CA row (fine, because WBITS even!)
	const size_t Q = I*q;     // half+1 bits in the CA
	double* const dps = malloc(Q*sizeof(double));
	TEST_ALLOC(dps);
	if (filtering) ca_dps(I,n,fca,dps,costab); else ca_dps(I,n,ca,dps,costab);
	scale(Q,dps,1.0/((double)m*(double)m));
	for (size_t i=0;i<Q;i+=q) dps[i] = NAN; // suppress S(0)
	const double dpsmax = max(Q,dps);
	const double dpsmin = min(Q,dps);
	FILE* const gp = gp_popen(NULL,"Discrete power spectrum",2480,1300);
	fprintf(gp,"# set size ratio -1\n");
	fprintf(gp,"unset xtics\n");
	fprintf(gp,"unset ytics\n");
	fprintf(gp,"set palette defined (%s)\n",gp_palette[0]);
	fprintf(gp,"set logs cb\n");
	fprintf(gp,"set cbr [1e-6:%g]\n",dpsmax);
	fprintf(gp,"set xr [+0.5:%g]\n",(double)q-0.5);
	fprintf(gp,"set yr [-0.5:%g]\n",(double)I-0.5);
	fprintf(gp,"plot '-' binary array=(%zu,%zu) flip=y with image not\n",q,I);
	gp_binary_write(gp,Q,dps,gpipw); // NOTE: if gpipw set, dps is now unusable!
	if (pclose(gp) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");
	printf("DPS max = %g, min = %g\n",dpsmax,dpsmin);
	free(dps);
}

void caana_automi
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const word_t* const fca,
	const int           filtering,
	const int           amice,
	const int           gpipw
)
{
	printf("calculating CA auto-MI ... "); fflush(stdout);
	const size_t m = n*WBITS; // bits in a CA row
	const size_t q = m/2+1;   // half+1 bits in a CA row (fine, because WBITS even!)
	const size_t Q = I*q;     // half+1 bits in the CA

	// calculate aut-MI
	double* const ami = malloc(Q*sizeof(double));
	TEST_ALLOC(ami);
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
	TEST_ALLOC(amik);
	double* const med  = malloc(q*sizeof(double));
	TEST_ALLOC(med);
	med[0] = NAN;
	for (size_t k=1;k<q;++k) {
		for (size_t r=0;r<I;++r) amik[r] = ami[q*r+k];
		med[k] = median(I,amik,NULL,0);
	}

	// display AMI heat map
	printf("AMI max = %g, min = %g\n",amimax,amimin);
	FILE* const gp = gp_popen(NULL,"Auto-MI: heat map",2480,1300);
	fprintf(gp,"# set size ratio -1\n");
	fprintf(gp,"unset xtics\n");
	fprintf(gp,"unset ytics\n");
	fprintf(gp,"set palette defined (%s)\n",gp_palette[0]);
	fprintf(gp,"set logs cb\n");
	fprintf(gp,"set cbr [1e-6:1]\n");
	fprintf(gp,"set xr [+0.5:%g]\n",(double)q-0.5);
	fprintf(gp,"set yr [-0.5:%g]\n",(double)I-0.5);
	fprintf(gp,"plot '-' binary array=(%zu,%zu) flip=y with image not\n",q,I);
	gp_binary_write(gp,Q,ami,gpipw); // NOTE: if gpipw set, ami is now unusable!
	if (pclose(gp) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");

	// display medians per distance box plot
	FILE* const gp1 = gp_popen(NULL,"Auto-MI: medians",2480,480);
	fprintf(gp1,"set xr [0.5:%g]\n",(double)q+0.5);
	fprintf(gp1,"set xtics 100\n");
	fprintf(gp1,"set xlabel \"distance\"\n");
	fprintf(gp1,"set ylabel \"median MI\"\n");
	fprintf(gp1,"set style fill solid 1\n");
	fprintf(gp1,"set grid\n");
	fprintf(gp1,"plot '-' u 1:2 w boxes not\n");
	for (size_t k=1;k<q;++k) fprintf(gp1,"%4zu  %12.8f f\n",k,med[k]);
	if (pclose(gp1) == EOF) PEEXIT("failed to close pipe to Gnuplot\n");

	free(med);
	free(amik);
	free(ami);
}
