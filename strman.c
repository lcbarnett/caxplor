#include <string.h>
#include <stdarg.h>

#include "strman.h"

sm_t sm_push(sm_t top)
{
	sm_t const oldtop = top;
	top = malloc(sizeof(sm_node_t));
	PASSERT(top != NULL,"memory allocation failed");
	top->prev = oldtop;
	top->str = NULL;
	return top;
}

sm_t sm_pop(sm_t top)
{
	sm_t const oldtop = top;
	top = oldtop->prev;
	free(oldtop->str);
	free(oldtop);
	return top;
}

char* sm_cpy(sm_t top, char* str)
{
	top = sm_push(top);
	top->str = malloc((strlen(str)+1)*sizeof(char));
	PASSERT(top->str != NULL,"memory allocation failed");
	strcpy(top->str,str);
	return top->str;
}

char* _sm_cat(sm_t top, char* str1, char* str2, ...)
{
	top = sm_push(top);
	va_list args;
	va_start(args,str2);
	size_t len = strlen(str1);
	for (char* str = str2; str!= NULL; str = va_arg(args,char*)) len += strlen(str);
	va_end(args);
	top->str = malloc((len+1)*sizeof(char));
	PASSERT(top->str != NULL,"memory allocation failed");
	strcpy(top->str,str1);
	va_start(args,str2);
	for (char* str = str2; str!= NULL; str = va_arg(args,char*)) strcat(top->str,str);
	va_end(args);
	return top->str;
}

void check_sm_printf_error(char* str)
{
	PASSERT(str != NULL,"memory allocation failed");
}
