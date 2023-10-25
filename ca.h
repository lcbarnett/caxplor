#ifndef CA_H
#define CA_H

#include "word.h"

/*********************************************************************/
/*                      CA (multi-word)                              */
/*********************************************************************/

void    ca_fprint      (const size_t I, const size_t n, const word_t* const ca, FILE* const fstream);
void    ca_print       (const size_t I, const size_t n, const word_t* const ca);
void    ca_fprints     (const size_t I, const size_t n, const word_t* const ca, FILE* const fstream);
void    ca_prints      (const size_t I, const size_t n, const word_t* const ca);
void    ca_part_count  (const size_t I, const size_t n, const word_t* const ca, const int P, const size_t nparts, partf_t* const ppw, const int sortem);
size_t  ca_period      (const size_t I, const size_t n, const word_t* const ca, const int B, const word_t* const rtab, int* const rot);

void    ca_rotl        (const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits);
void    ca_rotr        (const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int nbits);
void    ca_reverse     (const size_t I, const size_t n, word_t* const ca, const word_t* const caold);
void    ca_filter      (const size_t I, const size_t n, word_t* const ca, const word_t* const caold, const int B, const word_t* const rtab);

void    ca_run         (const size_t I, const size_t n, word_t* const ca, word_t* const cawrk, const int B, const word_t* const rtab, const int uto);

void    ca_dps         (const size_t I, const size_t n, const word_t* const ca, double* const dps, const double* const costab);
void    ca_autocov     (const size_t I, const size_t n, const word_t* const ca, double* const ac);
void    ca_automi      (const size_t I, const size_t n, const word_t* const ca, double* const ami);

#endif // CA_H
