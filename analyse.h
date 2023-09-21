#ifndef CAANALYSE
#define CAANALYSE

#include "word.h"

void caana_automi
(
	// sizes
	const size_t        n,
	const size_t        I,
	const size_t        q,
	const size_t        Q,
	// CAs
	const word_t* const ca,
	const word_t* const fca,
	// flags
	const int           filtering,
	const int           amice,
	// Gnuplot parameters
	const int           gpipw
);

#endif // CAANALYSE
