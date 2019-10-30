#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "barrier_types.h"
#include "utils.h"

#include "compress.h"

static int  ratio;
static char *alphabet;
static int  alphabet_size     = -1;
static int  orig_stringlength = -1;

void
ini_pack_em(barrier_options opt)
{
  /* 7:1 compression using base 2 encoding */
  /* 5:1 compression using base 3 encoding */
  /* 3:1 compression using base 4, 5 or 6 encoding */
  alphabet = (char *)space(strlen(opt.GRAPH) * sizeof(char));
  if (sscanf(opt.GRAPH, "Q%d,%s", &alphabet_size, alphabet) != 2) {
    fprintf(stderr, "error in opt.GRAPH\n");
    exit(777);
  }

  if (alphabet_size != strlen(alphabet))
    fprintf(stderr, "wrong alphabet size\n");

  orig_stringlength = strlen(opt.seq);

  switch (alphabet_size) {
    case 2:
      ratio = 7;
      break;
    case 3:
      ratio = 5;
      break;
    case 4:
      ratio = 3;
      break;
    case 5:
      ratio = 3;
      break;
    case 6:
      ratio = 3;
      break;
    default:
      fprintf(stderr, "Alphabet size neither 3, 4, 5 or 6\n");
      break;
  }
}


int
letter2num(char c)
{
  char *pos;

  /*  F  L  R  U  X  Y */
  /*  0  1  2  3  4  5 */
  if (c == '\0')
    return 0;

  pos = strchr(alphabet, c);
  if (pos == NULL)
    fprintf(stderr, "error in string\n");

  return (int)(pos - alphabet);
  /* attenetion: the first 3 chars of opt.GRAPH are Q<number>, */
  /* so do net use the characters 'Q', ',', and any number in your alphabet */
}


char *
pack_em(const char *string)
{
  int           i, j, l, pi;
  unsigned char *packed;

  l                 = strlen(string);
  orig_stringlength = l;
  packed            =
    (unsigned char *)calloc(1, ((l + ratio - 1) / ratio + 1) * sizeof(unsigned char));

  j = i = pi = 0;
  while (i < l) {
    register unsigned char p;
    for (p = pi = 0; pi < ratio; pi++) {
      p *= alphabet_size;  /* alphabet_size == base */
      p += letter2num(string[i]);
      if (i < l)
        i++;
    }
    packed[j++] = (unsigned char)p + 1;
  }
  packed[j] = '\0';
  return (char *)packed;
}


char *
unpack_em(const char *packed)
{
  int                 i, j, l;
  char                *struc;
  unsigned const char *pp;

  l     = strlen(packed);
  pp    = (unsigned char *)packed;
  struc = (char *)calloc(1, (l * ratio + 1) * sizeof(char));

  j = 0;
  for (i = j = 0; i < l; i++) {
    unsigned char p;
    int           k, c;

    p = pp[i] - 1;
    for (k = (ratio - 1); k >= 0; k--) {
      c             = p % alphabet_size;
      p             /= alphabet_size;
      struc[j + k]  = alphabet[c];
    }
    j += ratio;
  }

  while (j >= orig_stringlength)
    struc[j--] = '\0';

  return struc;
}
