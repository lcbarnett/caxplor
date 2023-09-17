#include "ca.h"
#include "utils.h"

/*********************************************************************/
/*                      CA (multi-word)                              */
/*********************************************************************/

void ca_fprint(const size_t I, const size_t n, const word_t* const ca, FILE* const fstream)
{
	for (const word_t* w=ca;w<ca+I*n;w+=n) {mw_fprint(n,w,fstream); fnewline(fstream);}
}

void ca_print(const size_t I, const size_t n, const word_t* const ca)
{
	ca_fprint(I,n,ca,stdout);
}

void ca_fprints(const size_t I, const size_t n, const word_t* const ca, FILE* const fstream)
{
	for (const word_t* w=ca;w<ca+I*n;w+=n) {mw_fprints(n,w,fstream); fnewline(fstream);}
}

void ca_prints(const size_t I, const size_t n, const word_t* const ca)
{
	ca_fprints(I,n,ca,stdout);
}

void ca_run(const size_t I, const size_t n, word_t* const ca, word_t* const cawrk, const int B, const word_t* const rtab, const int uto)
{
	for (word_t* w=ca+n;w<ca+I*n;w+=n) mw_filter(n,w,w-n,B,rtab);
	if (uto) {
		ASSERT(cawrk != NULL,"Need CA work buffer to unwrap!");
		mw_copy(I*n,cawrk,ca);
		ca_rotl(I,n,ca,cawrk,uto);
	}
}

void ca_filter(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int B, const word_t* const rtab)
{
	word_t* wnew = ca;
	for (const word_t* w=caold;w<caold+I*n;w+=n,wnew+=n) mw_filter(n,wnew,w,B,rtab);
}

size_t ca_period(const size_t I, const size_t n, const word_t* const ca, const int B, const word_t* const rtab, int* const rot)
{
	word_t mword1[n];
	word_t mword2[n];
	word_t* wold = mword1;
	word_t* wnew = mword2;
	mw_copy(n,wold,ca);
	for (size_t i=0;i<I;++i) {
		mw_filter(n,wnew,wold,B,rtab);
		*rot = mw_equiv(n,ca,wnew);
		if (*rot >= 0) return i+1;
		SWAP(word_t*,wnew,wold);
	}
	return I;
}

void ca_rotl(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits)
{
	size_t nb = (size_t)nbits;
	size_t b = 0;
	for (size_t r=0;r<I*n;r+=n,b+=nb) mw_rotl(n,ca+r,caold+r,b);
}

void ca_rotr(const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits)
{
	size_t nb = (size_t)nbits;
	size_t b = 0;
	for (size_t r=0; r<I*n; r+=n, b+=nb) mw_rotr(n,ca+r,caold+r,b);
}

void ca_reverse(const size_t I, const size_t n, word_t* const ca, const word_t* const caold)
{
	for (size_t r=0; r<I*n; r += n) mw_reverse(n,ca+r,caold+r);
}

void ca_part_count(const size_t I, const size_t n, const word_t* const ca, const int P, const size_t nparts, partf_t* const ppw, const int sortem)
{
	ASSERT(nparts == POW2(P),"number of particles must equal 2^(particle breadth)");
	word_t* const wcgrain = mw_alloc(n);
	const word_t* const caend = ca+I*n;
	for (word_t part=0;part<nparts;++part) {
		ulong ppwp = 0;
		for (const word_t* w=ca;w<caend;w+=n) {
			mw_zero(n,wcgrain);
			mw_parts(n,wcgrain,w,P,part);
			ppwp += (ulong)mw_nsetbits(n,wcgrain);
		}
		ppw[part].w = part;
		ppw[part].f = (double)ppwp/(double)I;
	}
	if (sortem) qsort(ppw,nparts,sizeof(partf_t),partf_comp);
	free(wcgrain);
}

void ca_dps(const size_t I, const size_t n, const word_t* const ca, double* const dps, const double* const costab)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	double* const dftre = calloc(q,sizeof(double));
	double* const dftim = calloc(q,sizeof(double));
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const dpsr = dps+q*row;
		mw_dft(n,car,dftre,dftim,dpsr,costab);
	}
	free(dftim);
	free(dftre);
}

void ca_autocov(const size_t I, const size_t n, const word_t* const ca, double* const ac)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const acr = ac+q*row;
		mw_autocov(n,car,acr);
	}
}

void ca_automi(const size_t I, const size_t n, const word_t* const ca, double* const ami)
{
	const size_t m = n*WBITS;
	const size_t q = m/2+1; // fine, because WBITS even!
	for (size_t row=0; row<I; ++row) {
		const word_t* const car = ca+n*row;
		double* const amir = ami+q*row;
		mw_automi(n,car,amir);
	}
}

/*********************************************************************/
/*                      bitmap stuff                                 */
/*********************************************************************/

gdImagePtr ca_image_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const int           foncol
)
{
	// NOTE: you should call gdImageDestroy(im) after using the image!

	const int imx = ppc*(int)n*WBITS; // pixels per row
	const int imy = ppc*(int)I;       // pixels per column

	// create a gd image object (with space for 1-pixel border)
	gdImagePtr const im = gdImageCreate(imx+2,imy+2);

	// cell on/off colours
	const int col0 = gdImageColorAllocate(im,255,255,255);          // off cell colour
	const int col1 = gdImageColorAllocate(im,0,0,foncol ? 127 : 0); // on  cell colour
	const int colb = gdImageColorAllocate(im,0,0,0);                // border colour

	// draw border, and set background to off colour
	gdImageRectangle(im,0,0,imx+1,imy+1,colb);   // 1-pixel border
	gdImageFilledRectangle(im,1,1,imx,imy,col0); // background

	// if no CA, return an empty image
	if (ca == NULL) return im;

	const word_t* w = ca;
	if (ppc == 1) {
		// draw a pixel at each on cell (reverse x direction, so left -> right is hi -> lo)
		for (int y=0; y<imy; ++y) { // for each row
			for (int x=0; x<imx; ++w) { // for each word
				for (word_t b=0; b<WBITS; ++b, ++x) { // for each bit in word
					if (BITON(*w,b)) gdImageSetPixel(im,imx-x,y+1,col1);
				}
			}
		}
	}
	else {
		// draw a rectangle at each on cell (reverse x direction, so left -> right is hi -> lo)
		const int revx = imx-ppc+1;
		for (int y=0; y<imy; y += ppc) { // for each row
			for (int x=0; x<imx; ++w) { // for each word
				for (word_t b=0; b<WBITS; ++b, x += ppc) { // for each bit in word
					if (BITON(*w,b)) gdImageFilledRectangle(im,revx-x,y+1,imx-x,y+ppc,col1);
				}
			}
		}
	}

	return im;
}

void ca_image_write(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const char*   const imfmt,
	FILE* const         imfp,
	const int           foncol
)
{
	// pointer to the gd image
	const gdImagePtr im = ca_image_create(I,n,ca,ppc,foncol);

	// write gd image to file in specified format
	if      (strcmp(imfmt,"png")==0) gdImagePng(im,imfp);
	else if (strcmp(imfmt,"bmp")==0) gdImageBmp(im,imfp,0);
	else if (strcmp(imfmt,"gif")==0) gdImageGif(im,imfp);
	else if (strcmp(imfmt,"jpg")==0||strcmp(imfmt,"jpeg")==0) gdImageJpeg(im,imfp,-1);
	else EEXIT("unhandled image format '%s'\n",imfmt);

	// free gd image memory
	gdImageDestroy(im); // free image buffer
}

void ca_image_write_file(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const char*   const imfmt,
	const char*   const imfile,
	const int           foncol
)
{
	// open the image file
	FILE* const imfp = fopen(imfile,"wb");
	if (imfp == NULL) PEEXIT("failed to open image file\n");

	// write the CA image
	ca_image_write(I,n,ca,ppc,imfmt,imfp,foncol);

	// close the image file
	if (fclose(imfp) == -1) PEEXIT(" failed to close image file\n");
}

XRectangle* ca_xrects_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const int           imx,
	const int           imy,
	const int           gap,
	const int           iright,
	int* const          nrects
)
{
	// NOTE: you should call free(rects) when done with them!

	*nrects = mw_nsetbits(I*n,ca); // number of "on" cells

	XRectangle* const rects = calloc(sizeof(XRectangle),(size_t)*nrects);
	PASSERT(rects != NULL,"memory allocation failed");

	const ushort uscpx = (ushort)ppc;
	const int    xroff = (iright ? 2*(imx+gap+2)-ppc : imx+gap+2-ppc);
	const int    yoff  = gap+2;

	// draw a rectangle at each on cell (reverse x direction, so left -> right is hi -> lo)
	const word_t* w = ca;
	size_t r = 0;
	for (int y=0; y<imy; y += ppc) { // for each row
		for (int x=0; x<imx;++w) { // for each word
			for (word_t b=0; b<WBITS; ++b, x += ppc) { // for each bit in word
				if (BITON(*w,b)) {
					rects[r].x      = (short)(xroff-x);
					rects[r].y      = (short)(y+yoff);
					rects[r].width  = uscpx;
					rects[r].height = uscpx;
					++r;
				}
			}
		}
	}

	return rects;
}

int ca_rects_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	XRectangle*   const rects,
	const int           ppc,
	const int           imx,
	const int           imy
)
{
	// Create a rectangle for each cell, on in first nrects, off in remainder.
	// NOTE: we reverse x direction, so left -> right is hi -> lo

	const int nrects = mw_nsetbits(I*n,ca); // number of "on" cells

	const ushort uscpx = (ushort)ppc;
	const int    xroff = imx+1-ppc;
	const int    yoff  = 1;

	const word_t* w = ca;
	int r1 = 0;
	int r2 = nrects;
	for (int y=0; y<imy; y += ppc) { // for each row
		for (int x=0; x<imx;++w) { // for each word
			for (word_t b=0; b<WBITS; ++b, x += ppc) { // for each bit in word
				if (BITON(*w,b)) {
					rects[r1].x      = (short)(xroff-x);
					rects[r1].y      = (short)(y+yoff);
					rects[r1].width  = uscpx;
					rects[r1].height = uscpx;
					++r1;
				}
				else {
					rects[r2].x      = (short)(xroff-x);
					rects[r2].y      = (short)(y+yoff);
					rects[r2].width  = uscpx;
					rects[r2].height = uscpx;
					++r2;
				}
			}
		}
	}

	return nrects;
}

#ifdef UNSAFE_ZPIXMAP

void ca_zpixmap_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	char* const         imdata,
	const int           ppc,
	const int           imx,
	const int           imy,
	const int           foncol
)
{
	// Build ZPixmap data for CA (risky, but neater and more efficient version)
	//
	// NOTE 1: this routine is pure evil and almost certainly not very portable,
	//         but it _is_ pretty damn efficient.
	//
	// NOTE 2: we DON'T reverse x direction, so left -> right is lo -> hi;
	//         if this bothers you, call ca_reverse() on the CA first.

	const uint32_t oncol  = foncol ? 0x000077 : 0x000000; // black, or dark blue for filter
	const uint32_t offcol = 0xFFFFFF;                     // white

	const word_t* w = ca;
	uint32_t* p = (uint32_t*)imdata; // Danger! Danger! We process pixels as (1st 24 bits of) 32-bit words

	if (ppc == 1) { // special case - avoid a whole lot of looping

		for (int y=0; y<imy; ++y) {     // for each cell row
			for (int x=0; x<imx; ++w) { // for each word
				for (word_t b=0; b<WBITS; ++b, ++x) *p++ = (BITON(*w,b) ? oncol : offcol); // for each pixel in row, set colour
			}
		}

	}
	else {

		const int xjmpf = imx-ppc;     // jump forward after pixel row
		const int xjmpb = ppc*(imx-1); // jump back after last pixel in pixels row

		for (int y=0; y<imy; y += ppc) { // for each pixel row
			for (int x=0; x<imx; ++w) {  // for each word
				for (word_t b=0; b<WBITS; ++b, x += ppc) { // for each bit in word (--> cell in row)
					if (BITON(*w,b)) {
						for (int u=0; u<ppc; ++u) {                  // for each row of pixels in cell
							for (int v=0; v<ppc; ++v) *p++ = oncol;  // for each pixel in row set colour
							p += xjmpf;                              // jump to top-left corner of next cell
						}
					}
					else {
						for (int u=0; u<ppc; ++u) {                  // for each row of pixels in cell
							for (int v=0; v<ppc; ++v) *p++ = offcol; // for each pixel in row set colour
							p += xjmpf;                              // jump to top-left corner of next cell
						}
					}
					p -= (x == xjmpf ? xjmpf : xjmpb); // jump to top-left corner of first cell in next cell row
				}
			}
		}

	}
}

#else // use safe(-ish) version

void ca_zpixmap_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	char* const         imdata,
	const int           ppc,
	const int           imx,
	const int           imy,
	const int           foncol
)
{
	// Build ZPixmap data for CA
	//
	// NOTE 1: this routine is still evil and probably not very portable either,
	//         but it's safer than the above, and still pretty efficient.
	//
	// NOTE 2: we DON'T reverse x direction, so left -> right is lo -> hi;
	//         if this bothers you, call ca_reverse() on the CA first.

	const char c = (char)255;
	const char d = (char)(foncol ? 127 : 0);

	const word_t* w = ca;
	char* p = imdata;

	if (ppc == 1) { // special case - avoid a whole lot of looping

		for (int y=0; y<imy; ++y) {     // for each cell row
			for (int x=0; x<imx; ++w) { // for each word
				for (word_t b=0; b<WBITS; ++b, ++x) { //  for each pixel in row, set colour
					if (BITON(*w,b)) {
						*p++ = d; // B
						*p++ = 0; // G
						*p++ = 0; // R
						*p++ = 0; // pad (probably unnecessary)
					}
					else {
						*p++ = c; // B
						*p++ = c; // G
						*p++ = c; // R
						*p++ = 0; // pad (probably unnecessary)
					}
				}
			}
		}

	}
	else {

		const int xjmpf = imx-ppc;     // x jump forward after cell row
		const int xjmpb = ppc*(imx-1); // x jump back after last cell in cell row
		const int pjmpf = 4*xjmpf;     // p jump forward after cell row
		const int pjmpb = 4*xjmpb;     // p jump back after last cell in cell row

		for (int y=0; y<imy; y += ppc) { // for each cell row
			for (int x=0; x<imx; ++w) {  // for each word
				for (word_t b=0; b<WBITS; ++b, x += ppc) { // for each bit in word (cell in row)
					if (BITON(*w,b)) {
						for (int u=0; u<ppc; ++u) {     // for each row of pixels in cell
							for (int v=0; v<ppc; ++v) { // for each pixel in row set colour
								*p++ = d; // B
								*p++ = 0; // G
								*p++ = 0; // R
								*p++ = 0; // pad (probably unnecessary)
							}
							p += pjmpf; // jump to top-left corner of next cell
						}
					}
					else {
						for (int u=0; u<ppc; ++u) {     // for each row of pixels in cell
							for (int v=0; v<ppc; ++v) { // for each pixel in row set colour
								*p++ = c; // B
								*p++ = c; // G
								*p++ = c; // R
								*p++ = 0; // pad (probably unnecessary)
							}
							p += pjmpf; // jump to top-left corner of next cell
						}
					}
					p -= (x == xjmpf ? pjmpf : pjmpb); // jump to top-left corner of first cell in next cell row
				}
			}
		}

	}
}

#endif // UNSAFE_ZPIXMAP
