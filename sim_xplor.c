#include <X11/Xlib.h>
#include <X11/Xutil.h>

#include "clap.h"
#include "screen_metrics.h"
#include "ca.h"
#ifdef HAVE_GD
	#include "cagd.h"
#endif
#include "caX11.h"
#include "rtab.h"
#include "strman.h"
#include "analyse.h"

void print_id(const rtl_t* const rule, const int filtering);

// Main "CA Explorer" simulation

int sim_xplor(int argc, char* argv[], int info)
{
	// CLAP (command-line argument parser). Default values
	// may be overriden on the command line as switches.
	//
	// Arg:   name     type     default       description
	puts("\n---------------------------------------------------------------------------------------");
	CLAP_CARG(nwords,  size_t,  0,            "number of words (or 0 for automatic)");
	CLAP_VARG(nrows,   size_t,  0,            "number of rows (or 0 for automatic)");
	CLAP_VARG(rsiz,    int,     5,            "CA rule size");
	CLAP_VARG(rlam,    double,  0.6,          "CA rule lambda");
	CLAP_CARG(rseed,   ulong,   0,            "CA rule random seed (or 0 for unpredictable)");
	CLAP_VARG(fsiz,    int,     0,            "filter rule size (or 0 for same as rule size)");
	CLAP_VARG(flam,    double,  0.8,          "filter rule lambda");
	CLAP_CARG(fseed,   ulong,   0,            "filter rule random seed (0 for unpredictable)");
	CLAP_CARG(iseed,   ulong,   0,            "initialisation random seed (0 for unpredictable)");
	CLAP_CARG(untwist, int,     1,            "untwist?");
	CLAP_CARG(irtfile, cstr,   "",            "input rtids file (empty to start with random rtid)");
	CLAP_CARG(ortfile, cstr,   "saved.rt",    "saved rtids file name");
	CLAP_CARG(utoff,   int,     0,            "untwist offset (or 0 for half-size)");
	CLAP_CARG(prff,    size_t,  1000000,      "period fast-forwrd");
	CLAP_CARG(pmax,    size_t,  100000,       "maximum period");
	CLAP_CARG(emmax,   int,     20,           "maximum sequence length for entropy calculation");
	CLAP_CARG(eiff,    int,     1,            "advance before entropy");
	CLAP_CARG(tmmax,   int,     14,           "maximum sequence length for DD calculation");
	CLAP_CARG(tiff,    int,     0,            "advance before DD calculation");
	CLAP_CARG(tlag,    int,     1,            "lag for DD calculation");
	CLAP_CARG(amice,   int,     0,            "auto-conditional entropy rather than auto-MI?");
	CLAP_CARG(ppc,     int,     1,            "cell display size in pixels");
	CLAP_CARG(gpx,     int,     32,           "horizontal gap in pixels");
	CLAP_CARG(gpy,     int,     32,           "vertical gap in pixels");
	CLAP_CARG(wdec,    int,     0,            "window decorations?");
	CLAP_CARG(gpdir,   cstr,   "/tmp",        "Gnuplot file directory");
	CLAP_CARG(gpipw,   int,     1,            "In-place binary write for Gnuplot");
#ifdef HAVE_GD
	CLAP_CARG(imdir,   cstr,   "/tmp",        "image file directory");
	CLAP_CARG(imfmt,   cstr,   "png",         "image format (png, bmp, gif, jpg/jpeg)");
#endif
	puts("---------------------------------------------------------------------------------------\n");

	if (info) return EXIT_SUCCESS; // display switches and return

	// get number of CA rows/cols/words to fit screen

	size_t nr, ncols, nrwords;
	get_ca_dims(ppc,gpx,gpy,&nr,&ncols,&nrwords,1);

	const size_t n = (nwords == 0 ? nrwords : nwords);
	const size_t I = (nrows  == 0 ? nr      : nrows);

	puts("\n---------------------------------------------------------------------------------------\n");

	fsiz = fsiz == 0 ? rsiz : fsiz;

	// pseudo-random number generators
	mt_t rrng, irng, frng;
	mt_seed(&rrng,rseed);
	mt_seed(&irng,iseed);
	mt_seed(&frng,fseed);

	const int uto = untwist ? utoff == 0 ? rsiz/2 : utoff : 0;

	// allocate CA storage
	const size_t ncawords = I*n;
	word_t* const ca  = mw_alloc(ncawords); // the CA
	word_t* const fca = mw_alloc(ncawords); // the filtered CA
	word_t* const wca = mw_alloc(ncawords); // working CA

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

	// the image data
	char* const imdata = im->data = malloc((uint)(im->bytes_per_line*imy));
	TEST_ALLOC(imdata);
	XInitImage(im);

	// control variables
	int filtering = 0; // xplorer mode
	int quit  = 0; // time to go

#ifdef HAVE_GD
	const size_t strbuflen = 100;
	char strbuf[strbuflen+1];
	int imseq = 0; // image sequence number
#endif

	// DFT tables
	double* const costab = dft_cstab_alloc(n*WBITS);

	const size_t mslen = 10;
	char modestr[] = "exploring";

	// terminal usagestr
	char usagestr[] =
		"Keys:\n"
		"-----\n\n"
		"h : display this help\n"
		"m : (or ESC) toggle CA/filter mode\n"
		"n : (or SPACE) new random CA/filter\n"
		"N : new CA/filter from user-supplied id\n"
		"d : (or DEL) delete CA/filter\n"
		"j : (or left-arrow) previous CA/filter\n"
		"k : (or right-arrow) next CA/filter\n"
		"J : first CA/filter\n"
		"K : last CA/filter\n"
		"c : change CA/filter size\n"
		"v : invert CA/filter\n"
		"f : forward CA one screen\n"
		"i : re-initialise CA\n"
		"E : calculate entropy of CA rule\n"
		"D : calculate dynamical dependence of CA/filter rules\n"
		"p : calculate CA period\n"
		"s : save CA/filter id to file\n"
#ifdef HAVE_GD
		"w : write CA image to file\n"
#endif
		"S : calculate CA spatial discrete power spectrum\n"
		"I : calculate CA spatial auto-MI\n"
		"q : exit program\n";
	printf("%s\n",usagestr);
	fflush(stdout);

	// initialise CA rule table list from file or random

	rtl_t* rule;
	if (irtfile[0] == '\0') { // no input rtid file
		rule = rtl_add(NULL,rsiz); // current CA rule: this should never be NULL!
		rt_randomise(rule->size,rule->tab,rlam,&rrng);
	}
	else { // have input rtid file
		printf("Reading rules and filters from '%s' ...\n",irtfile);
		FILE* const irtfs = fopen(irtfile,"r");
		if (irtfs == NULL) PEEXIT("failed to open input rtids file '%s'",irtfile);
		rule = rtl_fread(irtfs);
		ASSERT(rule != NULL,"No valid rtids found in input file!");
		if (fclose(irtfs) == -1) PEEXIT("failed to close input rtids file '%s'",irtfile);
		printf("Done\n\n");
	}
	printf("exploring : random CA : id = "); rt_print_id(rule->size,rule->tab);
	printf(", lambda = %6.4f\n",rt_lambda(rule->size,rule->tab));
	mw_randomise(n,ca,&irng);
	ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
	ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
	printf("%s : ",modestr);
	fflush(stdout);

	// saved CA rules file
	FILE* const ortfs = fopen(ortfile,"a");
	if (ortfs == NULL) PEEXIT("failed to open saved rtids file '%s'",ortfile);

	const int hlen = (emmax > tmmax ? emmax : tmmax)+1;
	double H [hlen];
	double Hf[hlen];
	double Tf[hlen];

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
//		if (key == 27) key = 'q'; // map ESC   to 'q' for quit
		if (key == 27) key = 'm'; // map ESC   to 'm' for toggle mode
		if (key == 32) key = 'n'; // map SPACE to 'n' for new
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
			if (filtering && rule->filt != NULL) {
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			strncpy(modestr,filtering ? "filtering" : "exploring",mslen);
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'n': // new random CA/filter

			if (filtering) {
				printf("random filter : ");
				rule->filt = rtl_add(rule->filt,fsiz);
				rt_randomise(rule->filt->size,rule->filt->tab,flam,&frng);
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				printf("random CA : ");
				rule = rtl_add(rule,rsiz);
				rt_randomise(rule->size,rule->tab,rlam,&rrng);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'N': // new user-supplied CA/filter

			if (filtering) {
				printf("enter filter id: "); // prompt for ftid
				fflush(stdout);
				word_t* const ftab = rt_read_id(&fsiz);
				printf("filtering : ");
				fflush(stdout);
				if (fsiz == -1) {
					printf("input is wrong length\n");
					break;
				}
				if (fsiz == -2) {
					printf("input contains non-hex characters\n");
					break;
				}
				printf("filter rule size = %d\n",fsiz);
				rule->filt = rtl_add(rule->filt,fsiz);
				rt_copy(rule->size,rule->filt->tab,ftab);
				free(ftab);
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
				printf("filtering : ");
				fflush(stdout);
			}
			else {
				printf("enter CA id: "); // prompt for rtid
				fflush(stdout);
				word_t* const rtab = rt_read_id(&rsiz);
				printf("exploring : ");
				fflush(stdout);
				if (rsiz == -1) {
					printf("input is wrong length\n");
					break;
				}
				if (rsiz == -2) {
					printf("input contains non-hex characters\n");
					break;
				}
				printf("CA rule size = %d\n",rsiz);
				rule = rtl_add(rule,rsiz);
				rt_copy(rule->size,rule->tab,rtab);
				free(rtab);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
				printf("exploring : ");
				fflush(stdout);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'd': // delete CA or filter from list

			if (filtering) {
				if (rule->filt == NULL) {
					printf("no filter to delete!\n");
					break;
				}
				printf("deleting filter : ");
				rule->filt = rtl_del(rule->filt);
				if (rule->filt == NULL) {
					ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
				}
				else {
					ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
					ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
				}
			}
			else {
				if (rule->prev == NULL && rule->next == NULL) {
					printf("last CA rule - won't delete!\n");
					break;
				}
				printf("deleting CA : ");
				rule = rtl_del(rule);
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'j': // back to previous CA or filter in list

			if (filtering) {
				if (rule->filt == NULL) {
					printf("previous filter : no filter!\n");
					break;
				}
				if (rule->filt->prev == NULL) {
					printf("previous filter : can't go back!\n");
					break;
				}
				printf("previous filter : ");
				rule->filt = rule->filt->prev;
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				if (rule->prev == NULL) {
					printf("previous CA : can't go back!\n");
					break;
				}
				printf("previous CA : ");
				rule = rule->prev;
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'k': // forward to next CA or filter in list

			if (filtering) {
				if (rule->filt == NULL) {
					printf("next filter : no filter!\n");
					break;
				}
				if (rule->filt->next == NULL) {
					printf("next filter : can't go forward!\n");
					break;
				}
				printf("next filter : ");
				rule->filt = rule->filt->next;
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				if (rule->next == NULL) {
					printf("next CA : can't go forward!\n");
					break;
				}
				printf("next CA : ");
				rule = rule->next;
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'J': // back to first CA or filter in list

			if (filtering) {
				if (rule->filt == NULL) {
					printf("first filter : no filter!\n");
					break;
				}
				if (rule->filt->prev == NULL) {
					printf("first filter : already at first!\n");
					break;
				}
				printf("first filter : ");
				while (rule->filt->prev != NULL) rule->filt = rule->filt->prev; // go to beginning of list
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				if (rule->prev == NULL) {
					printf("first CA : already at first!\n");
					break;
				}
				printf("first CA : ");
				while (rule->prev != NULL) rule = rule->prev; // go to beginning of list
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'K': // forward to last CA or filter in list

			if (filtering) {
				if (rule->filt == NULL) {
					printf("last filter : no filter!\n");
					break;
				}
				if (rule->filt->next == NULL) {
					printf("last filter : already at last!\n");
					break;
				}
				printf("last filter : ");
				while (rule->filt->next != NULL) rule->filt = rule->filt->next; // go to end of list
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				if (rule->next == NULL) {
					printf("last CA : already at last!\n");
					break;
				}
				printf("last CA : ");
				while (rule->next != NULL) rule = rule->next; // go to end of list
				mw_randomise(n,ca,&irng);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'c': // change CA/filter size (will apply to new CA/filter)

			if (filtering) {
				printf("enter new filter rule size : "); // prompt for size
				fflush(stdout);
				int newsize;
				const int ret = scanf("%d",&newsize);
				printf("exploring : ");
				if (ret == 0 || newsize == 0) {
					printf("bad size\n");
					break;
				}
				fsiz = newsize;
				printf("new filter rule size = %d\n",fsiz);
			}
			else {
				printf("enter new CA rule size : "); // prompt for size
				fflush(stdout);
				int newsize;
				const int ret = scanf("%d",&newsize);
				printf("exploring : ");
				if (ret == 0 || newsize == 0) {
					printf("bad size\n");
					break;
				}
				rsiz = newsize;
				printf("new CA rule size = %d\n",rsiz);
			}
			// no need to redisplay image
			break;

		case 'v': // invert CA/filter rule

			if (filtering) {
				if (rule->filt == NULL) {
					printf("inverting filter : no filter!\n");
					break;
				}
				printf("inverting filter : ");
				rt_invert(rule->filt->size,rule->filt->tab);
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
				flam = 1.0-flam;
			}
			else {
				printf("inverting CA : ");
				rt_invert(rule->size,rule->tab);
				ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
				rlam = 1.0-rlam;
			}
			print_id(rule,filtering);
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'f': // run CA forward one screen

			printf("fast-forward CA\n");
			fflush(stdout);
			mw_copy(n,ca,ca+(I-1)*n);
			ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
			if (filtering && rule->filt != NULL) {
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 'i': // re-initialise and rerun
			printf("re-initialise CA\n");
			mw_randomise(n,ca,&irng);
			ca_run(I,n,ca,wca,rule->size,rule->tab,uto);
			if (filtering && rule->filt != NULL) {
				ca_filter(I,n,fca,ca,rule->filt->size,rule->filt->tab);
				ca_zpixmap_create(I,n,fca,imdata,ppc,imx,imy,filtering);
			}
			else {
				ca_zpixmap_create(I,n,ca,imdata,ppc,imx,imy,filtering);
			}
			XPutImage(dis,win,gc,im,0,0,1,1,uimx,uimy);
			break;

		case 's': // save CA/filter rule id to file

			printf("saving CA ");
			rt_fprint_id(rule->size,rule->tab,ortfs);
			if (filtering && rule->filt != NULL) {
				printf("and filter rule ids\n");
				fputc(' ',ortfs);
				rt_fprint_id(rule->filt->size,rule->filt->tab,ortfs);
			}
			else {
				printf("rule id\n");
			}
			fputc('\n',ortfs);
			fflush(ortfs); // write to file
			// no need to redisplay image
			break;

#ifdef HAVE_GD

		case 'w': // write CA/filtered CA image to file

			++imseq;
			snprintf(strbuf,strbuflen,"%s/ca_%03d.%s",imdir,imseq,imfmt);
			printf("writing CA image to %s",strbuf);
			ca_image_write_file(I,n,ca,ppc,imfmt,strbuf,0);
			if (filtering && rule->filt != NULL) {
				snprintf(strbuf,strbuflen,"%s/fca_%03d.%s",imdir,imseq,imfmt);
				printf(", filtered image to %s",strbuf);
				fflush(stdout);
				ca_image_write_file(I,n,fca,ppc,imfmt,strbuf,0);
			}
			putchar('\n');
			// no need to redisplay image
			break;
#endif

		case 'p': // calculate CA period

			caana_period(n,I,ca,rule,prff,pmax);
			break;

		case 'E': // calculate entropy of CA rule

			printf("calculating CA/filter entropy");
			const size_t Se = POW2(emmax);
			TEST_RAM(Se*sizeof(uint64_t));
			uint64_t* const bine = malloc(Se*sizeof(uint64_t));
			TEST_ALLOC(bine);
			for (int m=0; m<hlen; ++m) H[m] = NAN;
			for (int m=rule->size; m<=emmax; ++m) H[m] = rt_entro(rule->size,rule->tab,m,eiff,bine)/(double)m;
			if (filtering && rule->filt != NULL) {
				for (int m=0; m<hlen; ++m) Hf[m] = NAN;
				for (int m=rule->filt->size; m<=emmax; ++m) Hf[m] = rt_entro(rule->filt->size,rule->filt->tab,m,eiff,bine)/(double)m;
			}
			free(bine);
			char gpename[] = "caentro";
			FILE* const gped = gp_dopen(gpename,gpdir);
			if (filtering && rule->filt != NULL) {
				printf(" rule entropy = %8.6f, filter entropy = %8.6f\n",H[emmax],Hf[emmax]);
				for (int m=0; m<hlen; ++m) fprintf(gped,"%d\t%g\t%g\n",m,H[m],Hf[m]);
			}
			else {
				printf(" rule entropy = %8.6f\n",H[emmax]);
				for (int m=0; m<hlen; ++m) fprintf(gped,"%d\t%g\n",m,H[m]);
			}
			if (fclose(gped) == -1) PEEXIT("failed to close Gnuplot data file\n");
			FILE* const gpec = gp_fopen(gpename,gpdir,NULL,"CA rule entropy",0,0);
			fprintf(gpec,"datfile = \"%s.dat\"\n",gpename);
			fprintf(gpec,"set title \"{/:Bold CA entropy}\\n\\nrule "); rt_fprint_id(rule->size,rule->tab,gpec); fprintf(gpec," ({/Symbol l} = %g)",rt_lambda(rule->size,rule->tab));
			if (rule->filt != NULL) {
				fprintf(gpec,", filter "); rt_fprint_id(rule->filt->size,rule->filt->tab,gpec); fprintf(gpec," ({/Symbol l} = %g)\"\n",rt_lambda(rule->filt->size,rule->filt->tab));
			}
			else {
				fprintf(gpec,"\"\n");
			}
			fprintf(gpec,"set xlabel \"CA length (bits)\"\n");
			fprintf(gpec,"set ylabel \"normalised entropy\"\n");
			fprintf(gpec,"set key right bottom Left rev\n");
			fprintf(gpec,"set grid\n");
			fprintf(gpec,"set xr [1:%d]\n",emmax);
			fprintf(gpec,"set yr [0:1]\n");
			fprintf(gpec,"set ytics 0.1\n");
			if (filtering && rule->filt != NULL) {
				fprintf(gpec,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy'\n");
			}
			else {
				fprintf(gpec,"plot datfile u 1:2 w lines t 'rule entropy'\n");
			}
			if (fclose(gpec) == -1) PEEXIT("failed to close Gnuplot command file\n");
			gp_fplot(gpename,gpdir);
			break;

		case 'D': // calculate dynamical dependence of CA/filter rules

			if (!filtering) {
				printf("not in filtering mode!\n");
				break;
			}
			if (rule->filt == NULL) {
				printf("no filter!\n");
				break;
			}
			printf("calculating CA/filter dynamical dependence");
			const size_t St = POW2(emmax);
			TEST_RAM(St*sizeof(uint64_t));
			uint64_t* const bint = malloc(St*sizeof(uint64_t));
			TEST_ALLOC(bint);
			const size_t S2t = POW2(2*tmmax);
			TEST_RAM(S2t*sizeof(uint64_t));
			uint64_t* const bin2t = malloc(S2t*sizeof(uint64_t));
			TEST_ALLOC(bin2t);
			for (int m=0; m<hlen; ++m) H[m] = NAN;
			for (int m=rule->size; m<=emmax; ++m) H[m] = rt_entro(rule->size,rule->tab,m,eiff,bint)/(double)m;
			for (int m=0; m<hlen; ++m) Hf[m] = NAN;
			for (int m=rule->filt->size; m<=emmax; ++m) Hf[m] = rt_entro(rule->filt->size,rule->filt->tab,m,eiff,bint)/(double)m;
			const int mmin = rule->size > rule->filt->size ? rule->size : rule->filt->size;
			for (int m=0; m<hlen; ++m) Tf[m] = NAN;
			for (int m=mmin; m<=tmmax; ++m) Tf[m] = rt_dd(rule->size,rule->tab,rule->filt->size,rule->filt->tab,m,tiff,tlag,bint,bin2t)/(double)m;
			free(bin2t);
			free(bint);
			printf(" rule entropy = %8.6f, filter entropy = %8.6f, DD = %8.6f\n",H[emmax],Hf[emmax],Tf[tmmax]);
			char gptname[] = "cadd";
			FILE* const gptd = gp_dopen(gptname,gpdir);
			for (int m=0; m<hlen; ++m) fprintf(gptd,"%d\t%g\t%g\t%g\n",m,H[m],Hf[m],Tf[m]);
			if (fclose(gptd) == -1) PEEXIT("failed to close Gnuplot data file\n");
			FILE* const gptc = gp_fopen(gptname,gpdir,NULL,"CA rule 1-lag Dynamical Dependence",0,0);
			fprintf(gptc,"datfile = \"%s.dat\"\n",gptname);
			fprintf(gptc,"set title \"{/:Bold CA dynamical dependence}\\n\\nrule "); rt_fprint_id(rule->size,rule->tab,gptc); fprintf(gptc," ({/Symbol l} = %g)",rt_lambda(rule->size,rule->tab));
			fprintf(gptc,", filter "); rt_fprint_id(rule->filt->size,rule->filt->tab,gptc); fprintf(gptc," ({/Symbol l} = %g)\"\n",rt_lambda(rule->filt->size,rule->filt->tab));
			fprintf(gptc,"set xlabel \"CA length (bits)\"\n");
			fprintf(gptc,"set ylabel \"normalised entropy\"\n");
			fprintf(gptc,"set key right bottom Left rev\n");
			fprintf(gptc,"set grid\n");
			fprintf(gptc,"set xr [1:%d]\n",emmax);
			fprintf(gptc,"set yr [0:1]\n");
			fprintf(gptc,"set ytics 0.1\n");
			fprintf(gptc,"plot datfile u 1:2 w lines t 'rule entropy', datfile u 1:3 w lines t 'filter entropy', datfile u 1:4 w lines t 'rule/filter DD'\n");
			if (fclose(gptc) == -1) PEEXIT("failed to close Gnuplot command file\n");
			gp_fplot(gptname,gpdir);
			break;

		case 'S': // calculate CA spatial discrete power spectrum

			caana_dps(n,I,ca,fca,filtering,costab,gpipw);
			break;

		case 'I': // calculate CA spatial auto-MI

			caana_automi(n,I,ca,fca,filtering,amice,gpipw);
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

	if (fclose(ortfs) == -1) PEEXIT("failed to close saved rtids file '%s'",ortfile);

	free(costab);

	rtl_free(rule);

	free(wca);
	free(fca);
	free(ca);

	return EXIT_SUCCESS;
}

void print_id(const rtl_t* const rule, const int filtering)
{
	printf("CA id = ");
	rt_print_id(rule->size,rule->tab);
	printf(", lambda = %6.4f",rt_lambda(rule->size,rule->tab));
	if (rule->filt != NULL && filtering) {
		printf(" : filter id = ");
		rt_print_id(rule->filt->size,rule->filt->tab);
		printf(", lambda = %6.4f",rt_lambda(rule->filt->size,rule->filt->tab));
	}
	putchar('\n');
}
