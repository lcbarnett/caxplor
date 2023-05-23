#ifndef CA_H
#define CA_H

#include <gd.h>
#include <X11/Xlib.h>

#include "word.h"

/*********************************************************************/
/*              rule table (double-linked) list                      */
/*********************************************************************/

typedef struct rtl_node {
	word_t*          rtab;
	struct rtl_node* prev;
	struct rtl_node* next;
	struct rtl_node* filt; // pointer to filter list
} rtl_t;

rtl_t*  rtl_add  (rtl_t* curr, const int B); // insert after
rtl_t*  rtl_del  (rtl_t* curr);
void    rtl_free (rtl_t* curr);

/*********************************************************************/
/*                      rule table                                   */
/*********************************************************************/

static inline void rt_copy(const int B, word_t* const rtdest, const word_t* const rtsrc)
{
	memcpy(rtdest,rtsrc,POW2(B)*sizeof(word_t));
}

static inline size_t rt_nsetbits(const int B, const word_t* const rtab)
{
	size_t b = 0;
	for (size_t r=0;r<POW2(B);++r) if (rtab[r]) ++b;
	return b;
}

static inline double rt_lambda(const int B, const word_t* const rtab) // Langton's lambda
{
	return (double)rt_nsetbits(B,rtab)/(double)POW2(B);
}

static inline size_t rt_nwords(const int B)
{
	return B > 6 ? POW2(B-6) : 1;
}

static inline void rt_randomise(const int B, word_t* const rtab, const double lam, mt_t* const prng)
{
	for (size_t r=0;r<POW2(B);++r) rtab[r] = (mt_rand(prng) < lam ? WONE : WZERO);
}

static inline void rt_invert(const int B, word_t* const rtab)
{
	for (size_t r=0;r<POW2(B);++r) rtab[r] = WONE-rtab[r];
}

word_t* rt_alloc       (const int B);

void    rt_randomb     (const int B, word_t* const rtab, const size_t b, mt_t* const prng);
void    rt_from_mwords (const int B, word_t* const rtab, const size_t nrtwords, const word_t* const rtwords);
int     rt_fread_id    (const int B, word_t* const rtab, FILE* const fstream);
int     rt_read_id     (const int B, word_t* const rtab);

size_t  rt_uwords      (const int B, const word_t* const rtab, const int m);
void    rt_to_mwords   (const int B, const word_t* const rtab, const size_t nrtwords, word_t* const rtwords);
void    rt_fprint      (const int B, const word_t* const rtab, FILE* const fstream);
void    rt_print       (const int B, const word_t* const rtab);
void    rt_fprint_id   (const int B, const word_t* const rtab, FILE* const fstream);
void    rt_print_id    (const int B, const word_t* const rtab);
void    rt_entro_hist  (const int B, const word_t* const rtab, const int m, const int iff, ulong* const bin);
double  rt_entro       (const int B, const word_t* const rtab, const int m, const int iff);
void    rt_trent1_hist (const int B, const word_t* const rtab, const int F, const word_t* const ftab, const int m, const int iff, const int ilag, ulong* const bin, ulong* const bin2);
double  rt_trent1      (const int B, const word_t* const rtab, const int F, const word_t* const ftab, const int m, const int iff, const int ilag);

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

void    ca_run         (const size_t I, const size_t n, word_t* const ca, const int B, const word_t* const rtab);

/*********************************************************************/
/*                      bitmap stuff                                 */
/*********************************************************************/

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

static const unsigned char hexchar[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

#endif // WORD_H
