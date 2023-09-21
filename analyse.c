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

void caana_entro
(
	const rtl_t* const rule,
	const int          filtering,
	const int          emmax,
	const int          eiff,
	const char*  const gpdir
)
{
		printf("calculating CA entropy");
		const int hlen = emmax+1;
		double H[hlen];
		for (int m=0; m<hlen; ++m) H[m] = NAN;
		double Hf[hlen];
		for (int m=0; m<hlen; ++m) Hf[m] = NAN;
		for (int m=rule->size; m<=emmax; ++m) {
			H[m] = rt_entro(rule->size,rule->tab,m,eiff)/(double)m;
			if (rule->filt != NULL) Hf[m] = rt_entro(rule->filt->size,rule->filt->tab,m,eiff)/(double)m;
			putchar('.'); fflush(stdout);
		}
		char gpename[] = "caentro";
		FILE* const gpd = gp_dopen(gpename,gpdir);
		if (filtering && rule->filt != NULL) {
			printf(" rule entropy = %8.6f, filter entropy = %8.6f\n",H[emmax],Hf[emmax]);
			for (int m=rule->size; m<hlen; ++m) fprintf(gpd,"%d\t%g\t%g\n",m,H[m],Hf[m]);
		}
		else {
			printf(" rule entropy = %8.6f\n",H[emmax]);
			for (int m=rule->size; m<hlen; ++m) fprintf(gpd,"%d\t%g\n",m,H[m]);
		}
		if (fclose(gpd) == -1) PEEXIT("failed to close Gnuplot data file\n");
		FILE* const gpc = gp_fopen(gpename,gpdir,NULL,"CA rule entropy",0,0);
		fprintf(gpc,"datfile = \"%s.dat\"\n",gpename);
		fprintf(gpc,"set title \"{/:Bold CA entropy}\\n\\nrule "); rt_fprint_id(rule->size,rule->tab,gpc); fprintf(gpc," ({/Symbol l} = %g)",rt_lambda(rule->size,rule->tab));
		if (rule->filt != NULL) {
			fprintf(gpc,", filter "); rt_fprint_id(rule->filt->size,rule->filt->tab,gpc); fprintf(gpc," ({/Symbol l} = %g)\"\n",rt_lambda(rule->filt->size,rule->filt->tab));
		}
		else {
			fprintf(gpc,"\"\n");
		}
		fprintf(gpc,"set xlabel \"CA length (bits)\"\n");
		fprintf(gpc,"set ylabel \"normalised entropy\"\n");
		fprintf(gpc,"set key right bottom Left rev\n");
		fprintf(gpc,"set grid\n");
		fprintf(gpc,"set xr [%d:%d]\n",rule->size,emmax);
		fprintf(gpc,"set yr [0:1]\n");
		fprintf(gpc,"set ytics 0.1\n");
		if (filtering && rule->filt != NULL) {
			fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy'\n");
		}
		else {
			fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy'\n");
		}
		if (fclose(gpc) == -1) PEEXIT("failed to close Gnuplot command file\n");
		gp_fplot(gpename,gpdir);
}

void caana_dd
(
	const rtl_t* const rule,
	const int          filtering,
	const int          emmax,
	const int          eiff,
	const int          tmmax,
	const int          tiff,
	const int          tlag,
	const char*  const gpdir
)
{
	printf("calculating CA/filter dynamical dependence");
	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	double H[hlen];
	for (int m=0; m<hlen; ++m) H[m] = NAN;
	double Hf[hlen];
	for (int m=0; m<hlen; ++m) Hf[m] = NAN;
	double Tf[hlen];
	for (int m=0; m<hlen; ++m) Tf[m] = NAN;
	for (int m=rule->size; m<=emmax; ++m) {
		H[m]  = rt_entro(rule->size,rule->tab,m,eiff)/(double)m;
		Hf[m] = rt_entro(rule->filt->size,rule->filt->tab,m,eiff)/(double)m;
		putchar('.'); fflush(stdout);
	}
	for (int m=rule->size; m<=tmmax; ++m) {
		Tf[m] = rt_trent1(rule->size,rule->tab,rule->filt->size,rule->filt->tab,m,tiff,tlag)/(double)m;
		putchar('.'); fflush(stdout);
	}
	printf(" rule entropy = %8.6f, filter entropy = %8.6f, DD = %8.6f\n",H[emmax],Hf[emmax],Tf[tmmax]);
	char gptname[] = "cadd";
	FILE* const gpd = gp_dopen(gptname,gpdir);
	for (int m=rule->size; m<hlen; ++m) fprintf(gpd,"%d\t%g\t%g\t%g\n",m,H[m],Hf[m],Tf[m]);
	if (fclose(gpd) == -1) PEEXIT("failed to close Gnuplot data file\n");
	FILE* const gpc = gp_fopen(gptname,gpdir,NULL,"CA rule 1-lag Dynamical Dependence",0,0);
	fprintf(gpc,"datfile = \"%s.dat\"\n",gptname);
	fprintf(gpc,"set title \"{/:Bold CA dynamical dependence}\\n\\nrule "); rt_fprint_id(rule->size,rule->tab,gpc); fprintf(gpc," ({/Symbol l} = %g)",rt_lambda(rule->size,rule->tab));
	fprintf(gpc,", filter "); rt_fprint_id(rule->filt->size,rule->filt->tab,gpc); fprintf(gpc," ({/Symbol l} = %g)\"\n",rt_lambda(rule->filt->size,rule->filt->tab));
	fprintf(gpc,"set xlabel \"CA length (bits)\"\n");
	fprintf(gpc,"set ylabel \"normalised entropy\"\n");
	fprintf(gpc,"set key right bottom Left rev\n");
	fprintf(gpc,"set grid\n");
	fprintf(gpc,"set xr [%d:%d]\n",rule->size,emmax);
	fprintf(gpc,"set yr [0:1]\n");
	fprintf(gpc,"set ytics 0.1\n");
	fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy', datfile u 1:4 w lines t 'rule/filter DD'\n");
	if (fclose(gpc) == -1) PEEXIT("failed to close Gnuplot command file\n");
	gp_fplot(gptname,gpdir);
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
