#include "rtab.h"
#include "utils.h"

/*********************************************************************/
/*              rule table (double-linked) list                         */
/*********************************************************************/

rtl_t* rtl_add(rtl_t* curr, const int size) // insert at end
{
	if (curr == NULL) { // empty list
		curr = malloc(sizeof(rtl_t));
		TEST_ALLOC(curr);
		curr->next = NULL;
		curr->prev = NULL;
	}
	else {
		while (curr->next != NULL) curr = curr->next; // go to end of list
		rtl_t* const oldcurr = curr;
		curr = malloc(sizeof(rtl_t));
		TEST_ALLOC(curr);
		curr->prev = oldcurr;
		curr->next = NULL;
		oldcurr->next = curr;
	}
	curr->filt = NULL;
	curr->size = size;
	curr->tab = rt_alloc(size);
	return curr;
}

rtl_t* rtl_del(rtl_t* curr)
{
	if (curr == NULL) return NULL; // nothing to delete
	rtl_t* const oldcurr = curr;
	if (oldcurr->next == NULL) { // end of list
		if (oldcurr->prev == NULL) { // only 1 item in list!
			curr = NULL;
		}
		else {
			curr = oldcurr->prev;
			curr->next = NULL;
		}
	}
	else {
		curr = oldcurr->next;
		if (oldcurr->prev == NULL) { // beginning of list
			curr->prev = NULL;
		}
		else {
			curr->prev = oldcurr->prev;
			oldcurr->prev->next = curr;
		}
	}
	free(oldcurr->tab);
	rtl_free(oldcurr->filt);
	return curr;
}

void rtl_free(rtl_t* curr)
{
	while (curr != NULL) curr = rtl_del(curr);
}

rtl_t* rtl_find(const rtl_t* const rule, const int size, const word_t* const tab)
{
	if (rule == NULL) return NULL;
	const rtl_t* r = rule;
	const size_t S = POW2(size);
	while (r->prev != NULL) r = r->prev; // go to beginning of list
	while (r != NULL) {
		if (size == r->size) {
			if (mw_equal(S,tab,r->tab)) return (rtl_t*)r;
		}
		r = r->next;
	}
	return NULL;
}

rtl_t* rtl_init(rtl_t* rule)
{
	if (rule == NULL) return NULL;
	while (rule->prev != NULL) rule = rule->prev; // go to beginning of rule list
	for (; rule->next != NULL; rule = rule->next) {
		if (rule->filt != NULL) { // have a filter list
			while (rule->filt->prev != NULL) rule->filt = rule->filt->prev; // go to beginning of filter list
		}
	}
	while (rule->prev != NULL) rule = rule->prev; // go to beginning of rule list again
	return rule;
}

int* rtl_nitems(const rtl_t* const rule, int* const nrules, int* const nfilts)
{
	// should call rtl_init(rule) first !!!
	*nrules = 0;
	if (rule == NULL) return NULL;
	for (const rtl_t* r = rule; r != NULL; r = r->next) ++(*nrules);
	int* const nfperr = calloc((size_t)*nrules,sizeof(int)); // zero-initialises
	TEST_ALLOC(nfperr);
	int* nfpr = nfperr;
	*nfilts = 0;
	for (const rtl_t* r = rule; r != NULL; r = r->next,++nfpr) {
		for (const rtl_t* f = r->filt; f != NULL; f = f->next)  {++(*nfpr); ++(*nfilts);}
	}
	return nfperr;
}

rtl_t* rtl_fread(FILE* rtfs)
{
	// A .rt (rtids) file should have lines of the form:
	//
	//     CA-rtid {filter-rtid}
	//
	// where {...} indicates optional. Blank lines and
	// whitespace in lines are ignored, and lines beginning
	// with # are taken as comments and ignored.

	rtl_t*  rule   = NULL;
	char*   line   = NULL;
	size_t  len    = 0;
	int     nlines = 0;
	ssize_t ilen;
	char    delimit[]=" \t\r\n\v\f"; // POSIX whitespace characters

	// read lines
	while ((ilen = getline(&line,&len,rtfs)) != -1) {

		++nlines;
		printf("%3d :",nlines);

		if (ilen > 0) line[ilen-1] = '\0'; // strip trailing newline

		const char* token = strtok(line,delimit); // first string: a rule id
		if (token == NULL) {
			printf(" empty line\n");
			continue;
		}
		if (token[0] == '#') { // comment!
			printf(" comment line\n");
			continue;
		}

		printf(" CA id = %s",token);
		int rsiz;
		word_t* const rtab = rt_sread_id(token,&rsiz);
		if (rsiz == -1) {
			printf(" - ERROR: bad size - skipped\n");
			continue;
		}
		if (rsiz == -2) {
			printf(" - ERROR: id contains non-hex characters - skipped\n");
			continue;
		}
		rtl_t* const rrule = rtl_find(rule,rsiz,rtab); // rule already read?
		if (rrule == NULL) {
			printf(" (new)");
			rule = rtl_add(rule,rsiz);
			rt_copy(rsiz,rule->tab,rtab);
		}
		else {
			printf(" (old)");
		}
		free(rtab);

		token = strtok(NULL,delimit); // second string: a filter id
		if (token == NULL) { // no filter id
			putchar('\n');
			continue;
		}

		printf(", filter id = %s",token);
		int fsiz;
		word_t* const ftab = rt_sread_id(token,&fsiz);
		if (fsiz == -1) {
			printf(" - ERROR: bad size - skipped\n");
			continue;
		}
		if (fsiz == -2) {
			printf(" - ERROR: id contains non-hex characters - skipped\n");
			continue;
		}
		rtl_t* const frule = rtl_find(rule->filt,fsiz,ftab); // rule filter already read?
		if (frule == NULL) {
			printf(" (new)");
			rule->filt = rtl_add(rule->filt,fsiz);
			rt_copy(fsiz,rule->filt->tab,ftab);
		}
		else {
			printf(" (old)");
		}
		free(ftab);

		token = strtok(NULL,delimit);
		if (token != NULL) { // extra token
			printf(" - WARNING: extra tokens ignored\n");
			continue;
		}

		putchar('\n');
    }
	fflush(stdout);

	free(line);

	rule = rtl_init(rule); // initialise rtl

	return rule;
}

/*********************************************************************/
/*                      rule table                                   */
/*********************************************************************/

word_t* rt_alloc(const int size)
{
	ASSERT(size<WBITS,"rule too wide");
	word_t* const tab = calloc(POW2(size),sizeof(word_t)); // note: zero initialises (this is fine)
	TEST_ALLOC(tab);
	return tab;
}

void rt_randomb(const int size, word_t* const tab, const size_t b, mt_t* const prng)
{
	const size_t S = POW2(size);
	ASSERT(b<=S,"too many bits");
	memset(tab,0,S*sizeof(word_t));
	for (size_t i=0; i<b; ++i) tab[i] = WONE;
	for (size_t i=0; i<b; ++i) {
		size_t j = i+RANDI(size_t,S-i,prng);
		const word_t tmp = tab[i];
		tab[i] = tab[j];
		tab[j] = tmp;
	}
}

size_t rt_uwords(const int size, const word_t* const tab,const int m)
{
	// Find the number of unique words generated by rule tab of breadth size
	// operating on all words of length m.
	const size_t M = POW2(m);
	const int WB  = WBITS-size;
	const int Wm  = WBITS-m;
	const int WBm = WB+m;
	if (m > WBITS) EEXIT("word too long!");
	word_t* const w = mw_alloc(M);
	for (word_t v=0;v<M;++v) {
		word_t wv = 0;
		int i = 0;
		for (;i<m-size+1;++i) PUTBIT(wv,i,tab[(v<<(WB-i))>>WB]);
		for (;i<m    ;++i) PUTBIT(wv,i,tab[((v<<Wm)>>(Wm+i))|((v<<(WBm-i))>>WB)]);
		w[v] = wv;
	}
	size_t n = 1;
	for (size_t v=1;v<M;++v) {
		size_t u = 0;
		for (;u<v;++u) if (w[u] == w[v]) break;
		if (u == v) ++n;
	}
	free(w);
	return n;
}

void rt_to_mwords(const int size, const word_t* const tab, const size_t nrtwords, word_t* const rtwords)
{
	ASSERT(nrtwords == rt_nwords(size),"Wrong number of words!");
	const size_t S = POW2(size);
	word_t* p = rtwords;
	*p = WZERO;
	int i = 0;
	for (size_t r=0;r<S;++r) {
		PUTBIT(*p,i,tab[r]);
		if (++i == WBITS) {*(++p) = WZERO; i = 0;}
	}
}

void rt_from_mwords(const int size, word_t* const tab, const size_t nrtwords, const word_t* const rtwords)
{
	ASSERT(nrtwords == rt_nwords(size),"Wrong number of words!");
	const size_t S = POW2(size);
	const word_t* p = rtwords;
	int i = 0;
	for (size_t r=0;r<S;++r) {
		tab[r] = BITON(*p,i);
		if (++i == WBITS) {++p; i = 0;}
	}
}

void rt_fprint(const int size, const word_t* const tab, FILE* const fstream)
{
	const size_t S = POW2(size);
	for (size_t r=0;r<S;++r) fprintf(fstream,"%3zu = "PBP08" -> %"PRIw"\n",r,PBI08(r),tab[r]);
}

void rt_print(const int size, const word_t* const tab)
{
	rt_fprint(size,tab,stdout);
}

void rt_fprint_id(const int size, const word_t* const tab, FILE* const fstream)
{
	const size_t S = POW2(size);
	word_t u = 0;
	int i = 0;
	for (size_t r=0;r<S;++r) {
		PUTBIT(u,i,tab[r]);
		if ((++i)%4 == 0) {
			fputc(hexchar[u],fstream);
			u = 0;
			i = 0;
		}
	}
}

void rt_print_id(const int size, const word_t* const tab)
{
	rt_fprint_id(size,tab,stdout);
}

size_t rt_sprint_id(const int size, const word_t* const tab, size_t sbuflen, char* const str)
{
	// Note: *concatenates* rtid to string; string buffer must be long enough to accommodate it!
	const size_t S = POW2(size);
	const size_t C = rt_hexchars(size);
	const size_t slen = strlen(str);
	ASSERT(slen+C < sbuflen,"string buffer too short!"); // note extra char for NUL terminator
	word_t u = 0;
	int i = 0;
	size_t c = slen;
	for (size_t r=0;r<S;++r) {
		PUTBIT(u,i,tab[r]);
		if ((++i)%4 == 0) {
			str[c++] = hexchar[u];
			u = 0;
			i = 0;
		}
	}
	str[c] = '\0'; // null terminate
	return C; // number of chars written not including NUL terminator
}

word_t* rt_fread_id(FILE* const fstream, int* const size)  // allocates rule table on sucess - remember to free!
{
	char* str = NULL;
	size_t len = 0;
	const ssize_t ilen = getline(&str,&len,fstream);
	ASSERT(ilen != -1,"Read failed.");
	str[ilen-1] = '\0'; // strip trailing newline
	word_t* const tab = rt_sread_id(str,size);
	free(str);
	return tab;
}

word_t* rt_read_id(int* const size) // allocates rule table on sucess - remember to free!
{
	return rt_fread_id(stdin,size);
}

word_t* rt_sread_id(const char* const str, int* const size) // allocates rule table on sucess - remember to free!
{
	const size_t len = strlen(str);
	*size = rt_hexsize(len);
	if (*size == -1) return NULL; // failure - ID is bad size
	size_t r = 0;
	word_t* const tab = rt_alloc(*size);
	const int nepc = *size == 1 ? 2 : 4; // number of table entries per char
	for (size_t c=0;c<len;++c) {
		const word_t u = hex2word(str[c]);
		if (u == 999) {
			free(tab);
			*size = -2; // failure - non-hex chars
			return NULL;
		}
		for (int i=0;i<nepc;++i) tab[r++] = BITON(u,i);
	}
	return tab; // success
}

double rt_entro( // Entropy for CA rule on sequence of length m after iff iterations

	const int           size,
	const word_t* const tab,
	const int           m,
	const int           iff,
	uint64_t*     const bin
)
{
	// Construct histogram

	const size_t S = POW2(m);
	for (size_t y=0; y<S; ++y) bin[y] = 0;
	for (word_t x=WZERO; x<S; ++x) {
		word_t y = x;
		for (int i=0; i<iff; ++i) y = wd_filter(m,y,size,tab); // advance CA (at least 1)
		++bin[y];
	}

	// Calculate entropy

	const double f = 1.0/(double)S;
	double* const p = (double* const)bin; // alias histogram as double array (!)
	for (size_t y=0; y<S; ++y) p[y] = f*(double)bin[y];
	const double H = entro2(S,p);
	return H;
}

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
)
{
	// Construct histograms

	const size_t S  = POW2(m);
	const size_t S2 = POW2(2*m);
	for (size_t y=0; y<S;  ++y) bin[y]  = 0;
	for (size_t y=0; y<S2; ++y) bin2[y] = 0;
	for (word_t x=WZERO; x<S; ++x) {
		word_t y = x;
		for (int i=0; i<iff; ++i) y = wd_filter(m,y,rsiz,rtab);  // advance CA (may be zero)
		const word_t u = wd_filter(m,y,fsiz,ftab);               // filter CA
		++bin[u];
		for (int i=0; i<ilag; ++i) y = wd_filter(m,y,rsiz,rtab); // advance CA (at least 1)
		const word_t v = wd_filter(m,y,fsiz,ftab);               // filter CA
		++bin2[u+S*v];
	}

	// Calculate entropies

	const double f = 1.0/(double)S;
	double* const p = (double* const)bin; // alias histogram as double array (!)
	for (size_t y=0; y<S; ++y) p[y] = f*(double)bin[y];
	const double H = entro2(S,p);
	double* const p2 = (double* const)bin2; // alias histogram as double array (!)
	for (size_t y2=0; y2<S2; ++y2) p2[y2] = f*(double)bin2[y2];
	const double H2 = entro2(S2,p2);
	return H2-H;
}
