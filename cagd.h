#ifndef CAGD_H
#define CAGD_H

#include <gd.h>

#include "word.h"

gdImagePtr ca_image_create(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const int           foncol
);

void ca_image_write(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const char*   const imfmt,
	FILE* const         imfp,
	const int           foncol
);

void ca_image_write_file(
	const size_t        I,
	const size_t        n,
	const word_t* const ca,
	const int           ppc,
	const char*   const imfmt,
	const char*   const imfile,
	const int           foncol
);

#endif // CAGD_H
