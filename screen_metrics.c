#include <X11/Xlib.h>

#include "word.h"
#include "utils.h"

void get_cpx(const double cmm, const int cgpx, int* const cpx, size_t* const I, const int verb)
{
	// return cell size in pixels and CA rows given cell size in mm

	Display* dis = XOpenDisplay(NULL);
	Screen*  scr = DefaultScreenOfDisplay(dis);

	const int xspx = scr->width;
	const int yspx = scr->height - cgpx;
	const int xsmm = scr->mwidth;
	const int ysmm = scr->mheight;

	const double xpxpmm = (double)xspx/(double)xsmm;
	const double ypxpmm = (double)yspx/(double)ysmm;

	const double cpxd = ypxpmm*cmm;
	const int    cpx1  = (int)round(cpxd);
	*cpx = cpx1 < 1 ? 1 : cpx1;
	if (cpx1 < 1) fprintf(stderr,"WARNING: cell size too small - increasing to 1 pixel = %.2fmm\n",1.0/ypxpmm);

	const double Id = (double)yspx/(double)*cpx;
	*I  = (size_t)Id;

	if (verb) {
		printf("Screen size px : %4d x %4d\n",xspx,yspx);
		printf("Screen size mm : %4d x %4d\n",xsmm,ysmm);
		printf("Pixels per  mm : %4.2f x %4.2f\n",xpxpmm,ypxpmm);
		printf("Cell px        : %g (%dis)\n",cpxd,*cpx);
		printf("rows           : %g (%lu)\n",Id,*I);
	}

	XCloseDisplay(dis);
}

void get_ca_dims(const int ppc, const int gpx, const int gpy, size_t* const nrows, size_t* const ncols, size_t* const nwords, const int verb)
{
	// return number of nrows to fit screen height, given pixels-per-cell ppc,
	// and a horizontal/vertical gap in pixels, gpx/gpy

	ASSERT(ppc >  0,"bad pixels-per-cell");
	ASSERT(gpx >= 0,"bad pixel x-gap");
	ASSERT(gpy >= 0,"bad pixel y-gap");

	Display* dis = XOpenDisplay(NULL);
	Screen*  scr = DefaultScreenOfDisplay(dis);

	const int xspx = scr->width  - gpx;
	const int yspx = scr->height - gpy;

	XCloseDisplay(dis);

	const double drows = (double)yspx/(double)ppc;
	const double dcols = (double)xspx/(double)ppc;
	*nrows = (size_t)drows;
	*ncols = (size_t)dcols;

	ASSERT(*nrows > 0,"oops... no rows!");
	ASSERT(*ncols > 0,"oops... no columns!");

	*nwords = (size_t)((double)*ncols/(double)WBITS);

	if (verb) {
		printf("Screen size : %4d x %4d\n",xspx,yspx);
		printf("Cell size   : %d\n",ppc);
		printf("Gap         : %d x %d\n",gpx,gpy);
		printf("cells       : %lu X %lu\n",*nrows,*ncols);
		printf("words       : %lu\n",*nwords);
	}
}

void SetWinNodeco(Display* const dis, Window win) // Voodoo
{
	typedef struct {
		unsigned long flags;
		unsigned long functions;
		unsigned long decorations;
		long input_mode;
		unsigned long status;
	} MwmHints;
	Atom mwmHintsProperty = XInternAtom(dis,"_MOTIF_WM_HINTS",0);
	MwmHints hints;
	hints.flags = 2; // (1L << 1)
	hints.decorations = 0;
	XChangeProperty(dis,win,mwmHintsProperty,mwmHintsProperty,32,PropModeReplace,(unsigned char*)&hints,5);
}
