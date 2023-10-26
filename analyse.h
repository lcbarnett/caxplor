#ifndef CAANALYSE
#define CAANALYSE

#include "word.h"
#include "rtab.h"

void caana_period
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const rtl_t*  const rule,
	const size_t        prff,
	const size_t        pmax
);

void caana_dps
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const word_t* const fca,
	const int           filtering,
	const double* const costab,
	const int           gpipw
);

void caana_automi
(
	const size_t        n,
	const size_t        I,
	const word_t* const ca,
	const word_t* const fca,
	const int           filtering,
	const int           amice,
	const int           gpipw
);

#endif // CAANALYSE
