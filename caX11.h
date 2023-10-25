#ifndef CAX11_H
#define CAX11_H

#include <X11/Xlib.h>

#include "word.h"

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
);

int ca_rects_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	XRectangle*   const rects,
	const int           ppc,
	const int           imx,
	const int           imy
);

void ca_zpixmap_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	char* const         imdata,
	const int           ppc,
	const int           imx,
	const int           imy,
	const int           foncol
);

#endif // CAX11_H
