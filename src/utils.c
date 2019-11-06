/*
 *           utils.c
 *
 *   c  Ivo L Hofacker and Walter Fontana
 *      Vienna RNA package
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <stdbool.h>
#include "config.h"
#ifdef WITH_DMALLOC
#include "dmalloc.h"
#endif

#include "utils.h"

unsigned short  xsubi[3];
int             cut_point = -1;

/*-------------------------------------------------------------------------*/

void *
space(size_t size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, size)) == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "SPACE: requested size: %ld\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    nrerror("SPACE allocation failure -> no memory");
  }

  return pointer;
}


#ifdef WITH_DMALLOC
#define space(S) calloc(1, (S))
#endif

/*------------------------------------------------------------------------*/

#ifndef WITH_DMALLOC
/* dmalloc.h #define's xrealloc */
void *
xrealloc(void   *p,
         size_t size)
{
  if (p == 0)
    return space(size);

  p = (void *)realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "xrealloc: requested size: %ld\n", size);
      nrerror("xrealloc allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    nrerror("xrealloc allocation failure -> no memory");
  }

  return p;
}


#endif
/*------------------------------------------------------------------------*/

void
nrerror(const char message[])                   /* output message upon error */


{
  fprintf(stderr, "\n%s\n", message);
  exit(-1);
}


/*------------------------------------------------------------------------*/
void
init_rand(void)
{
  time_t t;

  time(&t);
  xsubi[0]  = (unsigned short)t;
  xsubi[1]  = (unsigned short)(t >> 16);
  xsubi[2]  = 5246;
}


/*------------------------------------------------------------------------*/

double
urn(void)
/*
 * uniform random number generator; urn() is in [0,1]
 * uses a linear congruential library routine
 * 48 bit arithmetic
 */
{
  return erand48(xsubi);
}


/*------------------------------------------------------------------------*/

int
int_urn(int from,
        int to)
{
  return ((int)(urn() * (to - from + 1))) + from;
}


/*------------------------------------------------------------------------*/

void
filecopy(FILE *from,
         FILE *to)
{
  int c;

  while ((c = getc(from)) != EOF)
    putc(c, to);
}


/*-----------------------------------------------------------------*/

char *
time_stamp(void)
{
  time_t cal_time;

  cal_time = time(NULL);
  return ctime(&cal_time);
}


/*-----------------------------------------------------------------*/

char *
random_string(int         l,
              const char  symbols[])
{
  char  *r;
  int   i, rn, base;

  base  = strlen(symbols);
  r     = (char *)space(sizeof(char) * (l + 1));

  for (i = 0; i < l; i++) {
    rn    = (int)(urn() * base); /* [0, base-1] */
    r[i]  = symbols[rn];
  }
  r[l] = '\0';
  return r;
}


/*-----------------------------------------------------------------*/

int
hamming(const char  *s1,
        const char  *s2)
{
  int h = 0;

  for (; *s1 && *s2; s1++, s2++)
    if (*s1 != *s2)
      h++;

  return h;
}


/*-----------------------------------------------------------------*/

char *
get_line(FILE *fp)              /* reads lines of arbitrary length from fp */
{
  char s[512], *line, *cp;

  line = NULL;
  do {
    if (fgets(s, 512, fp) == NULL)
      break;

    cp = strchr(s, '\n');
    if (cp != NULL)
      *cp = '\0';

    if (line == NULL)
      line = space(strlen(s) + 1);
    else
      line = (char *)xrealloc(line, strlen(s) + strlen(line) + 1);

    strcat(line, s);
  } while (cp == NULL);

  return line;
}


char *
tokenize(char *line)
{
  char  *pos, *copy = NULL;
  int   cut = -1;

  if (line) {
    copy = (char *)space(strlen(line) + 1);
    (void)sscanf(line, "%s", copy);
    pos = strchr(copy, '&');
    if (pos) {
      cut = (int)(pos - copy) + 1;
      if (cut >= strlen(copy))
        cut = -1;

      if (strchr(pos + 1, '&'))
        nrerror("more than one cut-point in input");

      for (; *pos; pos++)
        *pos = *(pos + 1);               /* splice out the & */
    }

    if (cut > -1) {
      if (cut_point == -1) {
        cut_point = cut;
      } else if (cut_point != cut) {
        fprintf(stderr, "cut_point = %d cut = %d\n", cut_point, cut);
        nrerror("Two different cut points.");
      }
    }

    free(line);
  }

  return copy;
}


char *
costring(char *string)
{
  char  *ctmp;
  int   len;

  len   = strlen(string);
  ctmp  = (char *)space((len + 2) * sizeof(char));
  /* first sequence */
  (void)strncpy(ctmp, string, cut_point - 1);
  /* spacer */
  ctmp[cut_point - 1] = '&';
  /* second sequence */
  (void)strcat(ctmp, string + cut_point - 1);
  free(string);

  return ctmp;
}


static int
pack51(const char *struc,
       char       *packed,
       int        l)
{
  int i, j, pi;

  j = i = pi = 0;
  while (i < l) {
    register int p;
    for (p = pi = 0; pi < 5; pi++) {
      p *= 3;
      switch (struc[i]) {
        case '(':
        case '\0':
          break;
        case '.':
          p++;
          break;
        case ')':
          p += 2;
          break;
        default:
          nrerror("pack_structure: illegal character in structure");
      }
      if (i < l)
        i++;
    }
    packed[j++] = (unsigned char)(p + 1); /* never use 0, so we can use
                                           * strcmp()  etc. */
  }
  packed[j] = '\0';                       /* for str*() functions */
  return j;                               /* position of \0 */
}


char *
pack_structure(const char *struc)
{
  bool  isstar = false;
  int   key_endofstr, len;
  char  *key  = NULL;
  char  *s    = NULL;

  len = strlen(struc);
  key = (char *)space(((len + 4) / 5 + 2) * sizeof(char));
  s   = (char *)strdup(struc);
  if (s[len - 1] == '*') {
    /* we have a binding competent structure */
    isstar  = true;
    len     = len - 1;
    s[len]  = '\0';
  }

  key_endofstr = pack51(s, key, len);
  free(s);

  if (isstar) {
    key[key_endofstr]     = '\377';
    key[key_endofstr + 1] = '\0';
  }

  return key;
}


static int
unpack51(const char *packed,
         char       *struc,
         int        l)
{
  unsigned const char *pp;
  char                code[3] = {
    '(', '.', ')'
  };
  int                 i, j;

  pp = (const unsigned char *)packed;

  for (i = j = 0; i < l; i++) {
    register int p, c, k;

    p = (int)pp[i] - 1;
    for (k = 4; k >= 0; k--) {
      c             = p % 3;
      p             /= 3;
      struc[j + k]  = code[c];
    }
    j += 5;
  }
  struc[j--] = '\0';
  while (struc[j] == '(') /* strip trailing ( */
    struc[j--] = '\0';
  return j + 1;
}


/*-----------------------------------------------------------------*/
char *
unpack_structure(const char *packed)
{
  bool  isstar = false;
  int   struc_endofstr, len;
  char  *struc  = NULL;
  char  *p      = NULL;

  p   = (char *)strdup(packed);
  len = strlen(packed);
  if (p[len - 1] == '\377') {
    /* we have a binding competent structure */
    isstar  = true;
    len     = len - 1;
    p[len]  = '\0';
    struc   = (char *)space((len * 5 + 2) * sizeof(char));
  } else {
    struc = (char *)space((len * 5 + 2) * sizeof(char));
  }

  /* +2 above is 'cause we will eventually append a * in compute_rates() */
  struc_endofstr = unpack51(p, struc, len);
  free(p);
  if (isstar) {
    struc[struc_endofstr]     = '*';
    struc[struc_endofstr + 1] = '\0';
  }

  return struc;
}


/*
 * char *unpack_structure(const char *packed) {
 *   /\* 5:1 compression using base 3 encoding *\/
 *   int i,j,l;
 *   char *struc;
 *   unsigned const char *pp;
 *   char code[3] = {'(', '.', ')'};
 */

/*
 *   l = (int) strlen(packed);
 *   pp = (const unsigned char *) packed;
 *   struc = (char *) space((l*5+1)*sizeof(char));   /\* up to 4 byte extra *\/
 */

/*
 *   return struc;
 * }
 */

/*---------------------------------------------------------------------------*/

short *
make_pair_table(const char *structure)
{
  /* returns array representation of structure.
   * table[i] is 0 if unpaired or j if (i.j) pair.  */
  int   i, j, hx;
  int   length;
  short *stack;
  short *table;

  length    = strlen(structure);
  stack     = (short *)space(sizeof(short) * (length + 1));
  table     = (short *)space(sizeof(short) * (length + 2));
  table[0]  = length;

  for (hx = 0, i = 1; i <= length; i++) {
    switch (structure[i - 1]) {
      case '(':
        stack[hx++] = i;
        break;
      case ')':
        j = stack[--hx];
        if (hx < 0) {
          fprintf(stderr, "%s\n", structure);
          nrerror("unbalanced brackets in make_pair_table");
        }

        table[i]  = j;
        table[j]  = i;
        break;
      default:    /* unpaired base, usually '.' */
        table[i] = 0;
        break;
    }
  }
  if (hx != 0) {
    fprintf(stderr, "%s\n", structure);
    nrerror("unbalanced brackets in make_pair_table");
  }

  free(stack);
  return table;
}


/*---------------------------------------------------------------------------*/

int
bp_distance(const char  *str1,
            const char  *str2)
{
  /*
   * dist = {number of base pairs in one structure but not in the other}
   * same as edit distance with pair_open pair_close as move set
   */
  int   dist, i, l;
  short *t1, *t2;

  dist  = 0;
  t1    = make_pair_table(str1);
  t2    = make_pair_table(str2);

  l = (t1[0] < t2[0]) ? t1[0] : t2[0]; /* minimum of the two lengths */

  for (i = 1; i <= l; i++)
    if (t1[i] != t2[i]) {
      if (t1[i] > i)
        dist++;

      if (t2[i] > i)
        dist++;
    }

  free(t1);
  free(t2);
  return dist;
}
