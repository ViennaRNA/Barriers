/* Last changed Time-stamp: <2017-10-20 11:42:37 mtw> */
/* stapel.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

static char UNUSED  rcsid[] = "$Id: stapel.c,v 1.1 2001/04/05 08:00:57 ivo Exp $";
#define BASIS_SIZE 128

static char         **v       = NULL;
static int          len       = 0;
static int          stapelTop = 0;
static int          maxSize   = BASIS_SIZE;

/* oeffentliche funktionen */
void ini_stapel(int size);   /* changed pfs 03 2001 */


void push(char *form);


char *pop(void);


int get_top(void);


void reset_stapel(void);


void free_stapel(void);


/**/
void
ini_stapel(int size)
{
  int i;

  len = size + 3;
  v   = (char **)space(BASIS_SIZE * sizeof(char *));
  for (i = 0; i < BASIS_SIZE; i++)
    v[i] = (char *)space(len * sizeof(char));
  stapelTop = 0;
}


/**/
void
push(char *form)
{
  int i;

  if (stapelTop >= maxSize) {
    maxSize *= 2;
    v       = (char **)xrealloc(v, maxSize * sizeof(char *));
    for (i = stapelTop; i < maxSize; i++)
      v[i] = (char *)space(len * sizeof(char));
  }

  strcpy(v[stapelTop++], form);
}


/**/
void
reset_stapel(void)
{
  int i;

  if (stapelTop > BASIS_SIZE) {
    for (i = BASIS_SIZE; i < maxSize; i++)
      free(v[i]);
    v       = (char **)xrealloc(v, BASIS_SIZE * sizeof(char *));
    maxSize = BASIS_SIZE;
    for (i = 0; i < maxSize; i++)
      v[i][0] = '\0';
  }

  stapelTop = 0;
}


/**/
char *
pop(void)
{
  if (stapelTop == 0)
    return NULL;
  else
    return v[--stapelTop];
}


/**/
int
get_top(void)
{
  return stapelTop;
}


/**/
void
free_stapel(void)
{
  int i;

  for (i = 0; i < maxSize; i++)
    free(v[i]);
  free(v);
}
