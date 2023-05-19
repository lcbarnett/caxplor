#ifndef SM_H
#define SM_H

// A simple string manager which maintains a stack of
// C-string variables with printf-style formatting.
//
// Create a strman object with
//
//    sm_create(sm);
//
// You may then allocate formatted strings with, e.g.:
//
//    char* const mystr = sm_printf("x = %g, n = %d",x,n);
//
// as many times as required. The strman object must be
// destroyed with
//
//    sm_destroy(sm);
//
// to free all allocated strings and clean up the stack.

#include "utils.h"

// Internals

typedef struct sm_node {
	char* str;
	struct sm_node* prev;
} sm_node_t;

typedef sm_node_t* sm_t;

sm_t   sm_push (sm_t top);
sm_t   sm_pop  (sm_t top);
char*  sm_cpy  (sm_t top, char* str);
char* _sm_cat  (sm_t top, char* str1, char* str2, ...);

void check_sm_printf_error(char* str);

// User interface macros

#define sm_create(SM) \
	FILE* sm_devnull_fs = fopen("/dev/null","w"); \
	sm_t SM = NULL

#define sm_destroy(SM) \
	while (SM != NULL) SM = sm_pop(SM); \
	fclose(sm_devnull_fs)

#define sm_cat(top,str1, ...) _sm_cat(top,str1,__VA_ARGS__,NULL)

// The logic here is that we fprintf the formatted string to /dev/null to get the
// number of bytes to allocate for sprintf. Inefficient maybe, but it works :-)

#define sm_printf(top,fmt,...) \
	(top = sm_push(top), \
	top->str = malloc(((size_t)fprintf(sm_devnull_fs,fmt,__VA_ARGS__)+1)*sizeof(char)), \
	check_sm_printf_error(top->str), \
	sprintf(top->str,fmt,__VA_ARGS__), \
	top->str)

#endif // SM_H
