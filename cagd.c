#include "cagd.h"

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
