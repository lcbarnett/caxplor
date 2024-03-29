#ifndef RTAB_H
#define RTAB_H

#include "word.h"

/*********************************************************************/
/*              rule table (double-linked) list                      */
/*********************************************************************/

typedef struct rtl_node {
	int              size;
	word_t*          tab;
	struct rtl_node* prev;
	struct rtl_node* next;
	struct rtl_node* filt; // pointer to filter list
} rtl_t;

rtl_t*  rtl_add    (rtl_t* curr, const int size); // insert after
rtl_t*  rtl_del    (rtl_t* curr);
void    rtl_free   (rtl_t* curr);
rtl_t*  rtl_find   (const rtl_t* const rule, const int size, const word_t* const tab);
rtl_t*  rtl_init   (rtl_t* rule);
int*    rtl_nitems (const rtl_t* const rule, int* const nrules, int* const nfilts);
rtl_t*  rtl_fread  (FILE* rtfs);

/*********************************************************************/
/*                      rule table                                   */
/*********************************************************************/

static inline void rt_copy(const int size, word_t* const rtdest, const word_t* const rtsrc)
{
	memcpy(rtdest,rtsrc,POW2(size)*sizeof(word_t));
}

static inline size_t rt_nsetbits(const int size, const word_t* const tab)
{
	size_t b = 0;
	for (size_t r=0;r<POW2(size);++r) if (tab[r]) ++b;
	return b;
}

static inline double rt_lambda(const int size, const word_t* const tab) // Langton's lambda
{
	return (double)rt_nsetbits(size,tab)/(double)POW2(size);
}

static inline size_t rt_nwords(const int size)
{
	return size > 6 ? POW2(size-6) : 1;
}

static inline size_t rt_hexchars(const int size)
{
	return size > 2 ? POW2(size-2) : 1;
}

static inline int rt_hexsize(const size_t len)
{
	if (len == 0) return -1; // error
	if (len == 1) return 1;
	size_t size = 3, slen = 1;
	while ((slen <<= 1) < len) ++size;
	if (slen > len) return -1; // error
	return (int)size;
}

static inline void rt_randomise(const int size, word_t* const tab, const double lam, mt_t* const prng)
{
	for (size_t r=0;r<POW2(size);++r) tab[r] = (mt_rand(prng) < lam ? WONE : WZERO);
}

static inline void rt_invert(const int size, word_t* const tab)
{
	for (size_t r=0;r<POW2(size);++r) tab[r] = WONE-tab[r];
}

word_t* rt_alloc       (const int size);

void    rt_randomb     (const int size, word_t* const tab, const size_t b, mt_t* const prng);
void    rt_from_mwords (const int size, word_t* const tab, const size_t nrtwords, const word_t* const rtwords);
word_t* rt_fread_id    (FILE* const fstream, int* const size);   // allocates rule table on sucess - remember to free!
word_t* rt_read_id     (int* const size);                        // allocates rule table on sucess - remember to free!
word_t* rt_sread_id    (const char* const str, int* const size); // allocates rule table on sucess - remember to free!

size_t  rt_uwords      (const int size, const word_t* const tab, const int m);
void    rt_to_mwords   (const int size, const word_t* const tab, const size_t nrtwords, word_t* const rtwords);
void    rt_fprint      (const int size, const word_t* const tab, FILE* const fstream);
void    rt_fprint_id   (const int size, const word_t* const tab, FILE* const fstream);
void    rt_print       (const int size, const word_t* const tab);
size_t  rt_sprint_id   (const int size, const word_t* const tab, size_t sbuflen, char* const str);
void    rt_print_id    (const int size, const word_t* const tab);

double rt_entro( // Entropy for CA rule on sequence of length m after iff iterations
	const int           size,
	const word_t* const tab,
	const int           m,
	const int           iff,
	uint64_t*     const bin
);

double rt_dd( // dynamical dependence for CA/filter rules on sequence of length m after iff iterations, with lag ilag
	const int           rsiz,
	const word_t* const rtab,
	const int           fsiz,
	const word_t* const ftab,
	const int           m,
	const int           iff,
	const int           ilag,
	uint64_t*     const bin,
	uint64_t*     const bin2
);

static const char hexchar[] = {'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};

static inline word_t hex2word(const char c)
{
	return (word_t)((c >= '0') & (c <= '9') ? c-48 : (c >= 'A') & (c <= 'F') ? c-55 : 999);
}


#endif // RTAB_H
