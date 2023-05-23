#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "screen_metrics.h"
#include "ca.h"
#include "clap.h"
#include "strman.h"

void print_id(const int B, const rtl_t* const rtl, const int F, const int filtering)
{
	printf("CA id = ");
	rt_print_id(B,rtl->rtab);
	printf(", lambda = %6.4f",rt_lambda(B,rtl->rtab));
	if (rtl->filt != NULL && filtering) {
		printf(" : filter id = ");
		rt_print_id(F,rtl->filt->rtab);
		printf(", lambda = %6.4f",rt_lambda(F,rtl->filt->rtab));
	}
	putchar('\n');
}

int sim_xplor(int argc, char* argv[])
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_VARG(n,       size_t,  0,            "number of words (or 0 for automatic)");
	CLAP_CARG(B,       int,     5,            "CA rule breadth");
	CLAP_VARG(rlam,    double,  0.6,          "CA rule lambda");
	CLAP_CARG(rseed,   ulong,   0,            "CA rule random seed (or 0 for unpredictable)");
	CLAP_VARG(F,       int,     0,            "filter rule breadth (or 0 for same as rule breadth)");
	CLAP_VARG(flam,    double,  0.8,          "filter rule lambda");
	CLAP_CARG(fseed,   ulong,   0,            "filter rule random seed (0 for unpredictable)");
	CLAP_CARG(iseed,   ulong,   0,            "initialisation random seed (0 for unpredictable)");
	CLAP_CARG(untwist, int,     1,            "untwist?");
	CLAP_CARG(utoff,   int,     0,            "untwist offset (or 0 for half-B)");
	CLAP_CARG(prff,    size_t,  1000000,      "period fast-forwrd");
	CLAP_CARG(pmax,    size_t,  100000,       "maximum period");
	CLAP_CARG(emmax,   int,     20,           "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,    int,     1,            "advance before entropy");
	CLAP_CARG(tiff,    int,     0,            "advance before DD calculation");
	CLAP_CARG(tmmax,   int,     14,           "maximum sequence length for DD calculation");
	CLAP_CARG(tlag,    int,     1,            "lag for DD calculation");
	CLAP_CARG(ppc,     int,     1,            "cell display size in pixels");
	CLAP_CARG(gpx,     int,     32,           "horizontal gap in pixels");
	CLAP_CARG(gpy,     int,     32,           "vertical gap in pixels");
	CLAP_CARG(wdec,    int,     0,            "window decorations?");
	CLAP_CARG(gpdir,   cstr,   "/tmp",        "Gnuplot file directory");
	CLAP_CARG(imdir,   cstr,   "/tmp",        "image file directory");
	CLAP_CARG(imfmt,   cstr,   "png",         "image format (png, bmp, gif, jpg/jpeg)");
	puts("---------------------------------------------------------------------------------------\n");

	const size_t strbuflen = 100;
	char strbuf[strbuflen];

	// get number of CA rows/cols/words to fit screen

	size_t nrows, ncols, nrwords;
	get_ca_dims(ppc,gpx,gpy,&nrows,&ncols,&nrwords,1);

	n = (n == 0 ? nrwords : n);
	const size_t I = nrows;

	puts("\n---------------------------------------------------------------------------------------\n");

	F = F == 0 ? B : F;

	// pseudo-random number generators
	mt_t rrng, irng, frng;
	mt_seed(&rrng,rseed);
	mt_seed(&irng,iseed);
	mt_seed(&frng,fseed);

	// saved CA rules file
	char rtfile[] = "saved.rt";
	FILE* const rtfs = fopen(rtfile,"a");
	if (rtfs == NULL) PEEXIT("failed to open saved rtids file '%s'",rtfile);

	const int uto = untwist ? utoff == 0 ? B/2 : utoff : 0;

	// allocate CA storage
	const size_t nwords = I*n;
	word_t* const ca  = mw_alloc(nwords); // the CA
	word_t* const fca = mw_alloc(nwords); // the filtered CA
	word_t* const wca = mw_alloc(nwords); // working CA

	// window drawing constants
	const int  imx    = ppc*(int)n*WBITS; // pixels per row
	const int  imy    = ppc*(int)I;       // pixels per column
	const uint xwidth = (uint)(imx+2);
	const uint ywidth = (uint)(imy+2);
	const uint uimx   = (uint)imx;
	const uint uimy   = (uint)imy;

	// set up and map main window
	Display* const dis = XOpenDisplay(NULL);
	PASSERT(dis != NULL,"failed to open display");
	const int scr = DefaultScreen(dis);
	const ulong cblack = XBlackPixel(dis,scr);
	const ulong cwhite = XWhitePixel(dis,scr);
	const Window win = XCreateSimpleWindow(dis,DefaultRootWindow(dis),0,0,xwidth,ywidth,0,cblack,cwhite);
	const GC gc = XCreateGC(dis,win,0,NULL);
	if (!wdec) SetWinNodeco(dis,win);
	XStoreName(dis,win,"CA Xplorer");
	XSelectInput(dis,win,ExposureMask|KeyPressMask); // these are the events we listen for
	XMapWindow(dis, win);

	// set up image
	const uint depth = (uint)XDefaultDepth(dis,scr);
	XImage* const im = XCreateImage(dis,CopyFromParent,depth,ZPixmap,0,NULL,uimx,uimy,32,0);

	// the image data (solid red)
	im->data = malloc((uint)(im->bytes_per_line*imy));
	XInitImage(im);

	// control variables
	int filtering = 0; // xplorer mode
	int quit  = 0; // time to go
	int imseq = 0; // image sequence number

	// entropy storage
	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	double H[hlen];
	for (int m=0; m<hlen; ++m) H[m] = NAN;
	double Hf[hlen];
	for (int m=0; m<hlen; ++m) Hf[m] = NAN;
	double Tf[hlen];
	for (int m=0; m<hlen; ++m) Tf[m] = NAN;

	const size_t mslen = 10;
	char modestr[] = "exploring";

	// Gnuplot file streams
	FILE *gpd, *gpc;

	// terminal usagestr
	char usagestr[] =
		"Keys:\n"
		"-----\n\n"
		"h : display this help\n"
		"m : toggle CA/filter mode\n"
		"n : new random CA/filter\n"
		"N : new CA/filter from user-supplied id\n"
		"d : (or DEL) delete CA/filter\n"
		"j : (or left-arrow) previous CA/filter\n"
		"k : (or right-arrow) next CA/filter\n"
		"v : invert CA/filter\n"
		"f : forward CA one screen\n"
		"i : re-initialise CA\n"
		"e : plot entropy of CA rule\n"
		"t : plot 1-lag DD of CA rule and filter rule\n"
		"p : calculate CA period\n"
		"s : save CA/filter id to file\n"
		"w : write CA image to file\n"
		"q : (or ESC) exit program\n";
	printf("%s\n",usagestr);
	fflush(stdout);

	// initialise CA rule table list with random

	rtl_t* rtl = rtl_add(NULL,B); // this should never be NULL - there should laways be a CA rule!
	rt_randomise(B,rtl->rtab,rlam,&rrng);
	printf("exploring : random CA : id = "); rt_print_id(B,rtl->rtab);
	printf(", lambda = %6.4f\n",rt_lambda(B,rtl->rtab));
	mw_randomise(n,ca,&irng);
	ca_run(I,n,ca,wca,B,rtl->rtab,uto);
	ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
	printf("%s : ",modestr);
	fflush(stdout);

	// window event loop
	while (1) {

		// get next event in queue
		XEvent event;
		XNextEvent(dis,&event);

		// repaint on expose
		if (event.type == Expose) {
			// draw border
			XSetForeground(dis,gc,cblack);
			XDrawRectangle(dis,win,gc,0,0,xwidth-1,ywidth-1);
			// refresh CA
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			continue;
		}

		// we only take action on (some) key-presses
		if (event.type != KeyPress) continue;
		char key;
		KeySym keysym;
		int kret = XLookupString(&event.xkey,&key,1,&keysym,NULL);
		if (key == 27) key = 'q'; // map ESC to 'q' for quit
		switch (keysym) {
			case 65361: key = 'j'; kret = 1; break; // map left-arrow  to 'j' for prev
			case 65363: key = 'k'; kret = 1; break; // map right-arrow to 'k' for next
			case 65535: key = 'd'; kret = 1; break; // map DEL to 'd' for delete
		}
		if (kret != 1) continue;

		switch (key) {

		case 'm': // toggle CA/filter mode

			printf("switching mode : ");
			filtering = 1-filtering;
			if (filtering && rtl->filt != NULL) {
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			strncpy(modestr,filtering ? "filtering" : "exploring",mslen);
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'n': // new random CA/filter

			if (filtering) {
				printf("random filter : ");
				rtl->filt = rtl_add(rtl->filt,F);
				rt_randomise(F,rtl->filt->rtab,flam,&frng);
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				printf("random CA : ");
				rtl = rtl_add(rtl,B);
				rt_randomise(B,rtl->rtab,rlam,&rrng);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'N': // new user-supplied CA/filter

			if (filtering) {
				printf("enter filter id : "); // prompt for ftid
				fflush(stdout);
				word_t* const ftab = rt_alloc(F);
				const int res = rt_read_id(F,ftab);
				printf("filtering : ");
				fflush(stdout);
				if (res == 1) {
					printf("input is wrong length\n");
					free(ftab);
					break;
				}
				if (res == 2) {
					printf("input contains non-hex characters\n");
					free(ftab);
					break;
				}
				rtl->filt = rtl_add(rtl->filt,F);
				rt_copy(B,rtl->filt->rtab,ftab);
				free(ftab);
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				printf("enter CA id : "); // prompt for rtid
				fflush(stdout);
				word_t* const rtab = rt_alloc(B);
				const int res = rt_read_id(B,rtab);
				printf("exploring : ");
				fflush(stdout);
				if (res == 1) {
					printf("input is wrong length\n");
					free(rtab);
					break;
				}
				if (res == 2) {
					printf("input contains non-hex characters\n");
					free(rtab);
					break;
				}
				rtl = rtl_add(rtl,B);
				rt_copy(B,rtl->rtab,rtab);
				free(rtab);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'd': // delete CA or filter from list

			if (filtering) {
				if (rtl->filt == NULL) {
					printf("no filter to delete!\n");
					break;
				}
				printf("deleting filter : ");
				rtl->filt = rtl_del(rtl->filt);
				if (rtl->filt == NULL) {
					ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
				}
				else {
					ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
					ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
				}
			}
			else {
				if (rtl->prev == NULL && rtl->next == NULL) {
					printf("last CA rule - won't delete!\n");
					break;
				}
				printf("deleting CA : ");
				rtl = rtl_del(rtl);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'j': // back to previous CA or filter in list

			if (filtering) {
				if (rtl->filt == NULL) {
					printf("previous filter : no filter!\n");
					break;
				}
				if (rtl->filt->prev == NULL) {
					printf("previous filter : can't go back!\n");
					break;
				}
				printf("previous filter : ");
				rtl->filt = rtl->filt->prev;
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				if (rtl->prev == NULL) {
					printf("previous CA : can't go back!\n");
					break;
				}
				printf("previous CA : ");
				rtl = rtl->prev;
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'k': // forward to next CA or filter in list

			if (filtering) {
				if (rtl->filt == NULL) {
					printf("next filter : no filter!\n");
					break;
				}
				if (rtl->filt->next == NULL) {
					printf("next filter : can't go forward!\n");
					break;
				}
				printf("next filter : ");
				rtl->filt = rtl->filt->next;
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				if (rtl->next == NULL) {
					printf("next CA : can't go forward!\n");
					break;
				}
				printf("next CA : ");
				rtl = rtl->next;
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'v': // invert CA/filter rule

			if (filtering) {
				if (rtl->filt == NULL) {
					printf("inverting filter : no filter!\n");
					break;
				}
				printf("inverting filter : ");
				rt_invert(F,rtl->filt->rtab);
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
				flam = 1.0-flam;
			}
			else {
				printf("inverting CA : ");
				rt_invert(B,rtl->rtab);
				ca_run(I,n,ca,wca,B,rtl->rtab,uto);
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
				rlam = 1.0-rlam;
			}
			print_id(B,rtl,F,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'f': // run CA forward one screen

			printf("fast-forward CA\n");
			fflush(stdout);
			mw_copy(n,ca,ca+(I-1)*n);
			ca_run(I,n,ca,wca,B,rtl->rtab,uto);
			if (rtl->filt != NULL) {
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'i': // re-initialise and rerun

			printf("re-initialise CA\n");
			mw_randomise(n,ca,&irng);
			ca_run(I,n,ca,wca,B,rtl->rtab,uto);
			if (rtl->filt != NULL) {
				ca_filter(I,n,fca,ca,F,rtl->filt->rtab);
				ca_zpixmap_create(I,n,fca,im->data,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,im->data,ppc,imx,imy,filtering);
			}
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 's': // save CA/filter rule id to file

			printf("saving CA ");
			fprintf(rtfs,"B =%2d, id = ",B);
			rt_fprint_id(B,rtl->rtab,rtfs);
			if (rtl->filt != NULL) {
				printf("and filter rule ids\n");
				fflush(stdout);
				fprintf(rtfs," : F =%2d, id = ",F);
				rt_fprint_id(F,rtl->filt->rtab,rtfs);
			}
			else {
				printf("rule id\n");
			}
				fputc('\n',rtfs);
			break;

		case 'w': // write CA/filtered CA image to file

			++imseq;
			snprintf(strbuf,strbuflen,"%s/ca_%03d.%s",imdir,imseq,imfmt);
			printf("writing CA image to %s",strbuf);
			ca_image_write_file(I,n,ca,ppc,imfmt,strbuf,0);
			if (rtl->filt != NULL) {
				snprintf(strbuf,strbuflen,"%s/fca_%03d.%s",imdir,imseq,imfmt);
				printf(", filtered image to %s",strbuf);
				fflush(stdout);
				ca_image_write_file(I,n,fca,ppc,imfmt,strbuf,0);
			}
			putchar('\n');
			break;

		case 'p': // calculate CA period

			printf("calculating CA period... "); fflush(stdout);
			mw_copy(nwords,wca,ca); // copy to working CA
			mw_run(prff,n,wca,B,rtl->rtab);
			int prot;
			const size_t period = ca_period(pmax,n,wca,B,rtl->rtab,&prot);
			if (period == pmax) {
				printf("> %zu iterations\n",pmax);
			}
			else {
				printf("%zu iterations, twist = %d\n",period,prot);
			}
			break;

		case 'e': // calculate CA/filter entropy

			printf("calculating CA entropy");
			for (int m=B; m<=emmax; ++m) {
				H[m] = rt_entro(B,rtl->rtab,m,eiff)/(double)m;
				if (rtl->filt != NULL) Hf[m] = rt_entro(F,rtl->filt->rtab,m,eiff)/(double)m;
				putchar('.'); fflush(stdout);
			}
			char gpename[] = "caentro";
			gpd = gp_dopen(gpename,gpdir);
			if (rtl->filt != NULL) {
				printf(" rule entropy = %8.6f, filter entropy = %8.6f\n",H[emmax],Hf[emmax]);
				for (int m=B; m<hlen; ++m) fprintf(gpd,"%d\t%g\t%g\n",m,H[m],Hf[m]);
			}
			else {
				printf(" rule entropy = %8.6f\n",H[emmax]);
				for (int m=B; m<hlen; ++m) fprintf(gpd,"%d\t%g\n",m,H[m]);
			}
			if (fclose(gpd) == -1) PEEXIT("failed to close Gnuplot data file\n");
			gpc = gp_fopen(gpename,gpdir,NULL);
			fprintf(gpc,"datfile = \"%s.dat\"\n",gpename);
			fprintf(gpc,"set title \"{/:Bold CA entropy}\\n\\nrule "); rt_fprint_id(B,rtl->rtab,gpc); fprintf(gpc," ({/Symbol l} = %g)",rt_lambda(B,rtl->rtab));
			if (rtl->filt != NULL) {
				fprintf(gpc,", filter "); rt_fprint_id(F,rtl->filt->rtab,gpc); fprintf(gpc," ({/Symbol l} = %g)\"\n",rt_lambda(F,rtl->filt->rtab));
			}
			else {
				fprintf(gpc,"\"\n");
			}
			fprintf(gpc,"set xlabel \"CA length (bits)\"\n");
			fprintf(gpc,"set ylabel \"normalised entropy\"\n");
			fprintf(gpc,"set key right bottom Left rev\n");
			fprintf(gpc,"set grid\n");
			fprintf(gpc,"set xr [%d:%d]\n",B,emmax);
			fprintf(gpc,"set yr [0:1]\n");
			fprintf(gpc,"set ytics 0.1\n");
			if (rtl->filt != NULL) {
				fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy'\n");
			}
			else {
				fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy'\n");
			}
			if (fclose(gpc) == -1) PEEXIT("failed to close Gnuplot command file\n");
			gp_fplot(gpename,gpdir,NULL);
			break;

		case 't': // calculate CA/filter 1-lag DD

			if (rtl->filt == NULL) {
				printf("no filter!\n");
				break;
			}
			printf("calculating CA/filter dynamical dependence");
			for (int m=B; m<=emmax; ++m) {
				H[m]  = rt_entro(B,rtl->rtab,m,eiff)/(double)m;
				Hf[m] = rt_entro(F,rtl->filt->rtab,m,eiff)/(double)m;
				putchar('.'); fflush(stdout);
			}
			for (int m=B; m<=tmmax; ++m) {
				Tf[m] = rt_trent1(B,rtl->rtab,F,rtl->filt->rtab,m,tiff,tlag)/(double)m;
				putchar('.'); fflush(stdout);
			}
			printf(" rule entropy = %8.6f, filter entropy = %8.6f, DD = %8.6f\n",H[emmax],Hf[emmax],Tf[tmmax]);
			char gptname[] = "cadd";
			gpd = gp_dopen(gptname,gpdir);
			for (int m=B; m<hlen; ++m) fprintf(gpd,"%d\t%g\t%g\t%g\n",m,H[m],Hf[m],Tf[m]);
			if (fclose(gpd) == -1) PEEXIT("failed to close Gnuplot data file\n");
			gpc = gp_fopen(gptname,gpdir,NULL);
			fprintf(gpc,"datfile = \"%s.dat\"\n",gptname);
			fprintf(gpc,"set title \"{/:Bold CA dynamical dependence}\\n\\nrule "); rt_fprint_id(B,rtl->rtab,gpc); fprintf(gpc," ({/Symbol l} = %g)",rt_lambda(B,rtl->rtab));
			fprintf(gpc,", filter "); rt_fprint_id(F,rtl->filt->rtab,gpc); fprintf(gpc," ({/Symbol l} = %g)\"\n",rt_lambda(F,rtl->filt->rtab));
			fprintf(gpc,"set xlabel \"CA length (bits)\"\n");
			fprintf(gpc,"set ylabel \"normalised entropy\"\n");
			fprintf(gpc,"set key right bottom Left rev\n");
			fprintf(gpc,"set grid\n");
			fprintf(gpc,"set xr [%d:%d]\n",B,emmax);
			fprintf(gpc,"set yr [0:1]\n");
			fprintf(gpc,"set ytics 0.1\n");
			fprintf(gpc,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy', datfile u 1:4 w lines t 'rule/filter DD'\n");
			if (fclose(gpc) == -1) PEEXIT("failed to close Gnuplot command file\n");
			gp_fplot(gptname,gpdir,NULL);
			break;

		case 'h': // display usage

			printf("help\n\n%s\n",usagestr);
			break;

		case 'q': // quit

			printf("quit\n");
			quit = 1;
			break;

		default:
			printf("unrecognised action\n");

		}

		if (quit) break;

		printf("%s : ",modestr);
		fflush(stdout);
	}

	// clean up

	XDestroyImage(im);
	XDestroyWindow(dis,win);
	XCloseDisplay(dis);

	if (fclose(rtfs) == -1) PEEXIT("failed to close saved rtids file '%s'",rtfile);

	rtl_free(rtl);

	free(wca);
	free(fca);
	free(ca);

	return EXIT_SUCCESS;
}
