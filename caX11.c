#include "caX11.h"

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
