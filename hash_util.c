/* Last changed Time-stamp: <2017-10-02 11:58:00 mtw> */
/* hash_util.c */

#include "config.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "utils.h"
#include "hash_util.h"

#define   PRIVATE   static
#define   PUBLIC

static char UNUSED rcsid[] = "$Id: hash_util.c,v 1.9 2006/07/17 09:22:31 xtof Exp $";

/* modify hash_f(), hash_comp() and the typedef of hash_entry in hash_utils.h
 * to suit your application */

PUBLIC void *lookup_hash(void *x);


PUBLIC int write_hash(void *x);


PUBLIC void delete_hash(void *x);


PUBLIC void kill_hash();


PUBLIC void initialize_hash();


PUBLIC int hash_comp(void *x,
                     void *y);


inline PRIVATE unsigned hash_f(void *x);


/* HASHBITS usually defined via configure and config.h */
#ifndef HASHBITS
#define HASHBITS 24
#endif
#define HASHSIZE (((unsigned long)1 << HASHBITS) - 1)
/* on 64bit machines: unsigned long num = 1<<63; */

/* #define HASHSIZE 67108864 -1 */ /* 2^26 -1   must be power of 2 -1 */
/* #define HASHSIZE 33554432 -1 */ /* 2^25 -1   must be power of 2 -1 */
/* #define HASHSIZE 16777216 -1 */ /* 2^24 -1   must be power of 2 -1 */
/* #define HASHSIZE 4194304 -1  */ /* 2^22 -1   must be power of 2 -1 */

PRIVATE void          *hashtab[HASHSIZE + 1];

PUBLIC unsigned long  collisions = 0;

//static unsigned long  hashfillmax = HASHSIZE * 2. / 3.;
static unsigned long  hashfillmax = HASHSIZE;

unsigned long Number_of_hash_entries = 0;

/* ----------------------------------------------------------------- */

/* stolen from perl source */
char coeff[] = {
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1,
  61, 59, 53, 47, 43, 41, 37, 31, 29, 23, 17, 13, 11, 7, 3, 1
};

/* key must not be longer than 128 */
inline
PRIVATE unsigned
old_hash_f(void *x)
{
  register char *s;
  register int  i;
  register int  hash;

  s = ((hash_entry *)x)->structure;
  for (i = 0, hash = 0;
       /* while */ *s;
       s++, i++, hash *= 5)
    hash += *s * coeff[i];
  /* printf("%7d\t", hash); */
  return (hash) & (HASHSIZE);   /* modulo HASHSIZE for normalization */
}


PUBLIC int
hash_comp(void  *x,
          void  *y)
{
  return strcmp(((hash_entry *)x)->structure, ((hash_entry *)y)->structure);
}


/* ----------------------------------------------------------------- */

PUBLIC void *
lookup_hash(void *x)                 /* returns NULL unless x is in the hash */
{
  unsigned int hashval;

  hashval = hash_f(x);
  /* xtof poset debug ! */
#ifdef _DEBUG_HASH_
  fprintf(stderr,
          "lookup %s => %d\n",
          ((hash_entry *)x)->structure,
          hashval);
#endif
  if (hashtab[hashval] == NULL)
    return NULL;
  unsigned long hash_value_rejections = 0;
  while (hashtab[hashval]) {
    if (hash_comp(x, hashtab[hashval]) == 0)
      return hashtab[hashval];
    hash_value_rejections++;
    if (hash_value_rejections >= hashfillmax) {
      return NULL;
    }
    hashval = ((hashval + 1) & (HASHSIZE));
  }
  return NULL;
}


/* ----------------------------------------------------------------- */

PUBLIC int
write_hash(void *x)               /* returns 1 if x already was in the hash */
{
  unsigned int hashval;

  hashval = hash_f(x);
#ifdef _DEBUG_HASH_
  fprintf(stderr,
          "write  %s => %d\n",
          ((hash_entry *)x)->structure,
          hashval);
#endif
  unsigned long hash_value_rejections = 0;
  while (hashtab[hashval]) {
    if (hash_comp(x, hashtab[hashval]) == 0)
      return 1;

    hashval = ((hashval + 1) & (HASHSIZE));
    collisions++;
    hash_value_rejections++;
    if (hash_value_rejections >= hashfillmax) {
      /* die if # of hash entries exceeds 0.75*HASHSIE */
      float percentage = 100.0f * ((float)hashfillmax/(float)HASHSIZE);
      fprintf(stderr, "Warning: The hash is filled up to %.0f%% of hashsize (%ld | %ld)! Structure %s has been omitted!\n", percentage, hashfillmax, HASHSIZE, ((hash_entry *)x)->structure);
      return -1;
    }
  }
  hashtab[hashval] = x;
  Number_of_hash_entries++;
  return 0;
}


/* ----------------------------------------------------------------- */

PUBLIC void
initialize_hash()
{
}


/* ----------------------------------------------------------------- */

PUBLIC void
kill_hash()
{
  unsigned int i;

  for (i = 0; i < HASHSIZE + 1; i++) {
    if (hashtab[i]) {
      free(((hash_entry *)hashtab[i])->structure);
      /*free (hashtab[i]);*/
      hashtab[i] = NULL;
    }
  }
  Number_of_hash_entries = 0;
}


/* ----------------------------------------------------------------- */

PUBLIC void
delete_hash(void *x)               /* doesn't work in case of collsions */
{
  /* doesn't free anything ! */
  unsigned int hashval;

  unsigned long hash_value_rejections = 0;
  hashval = hash_f(x);
  while (hashtab[hashval]) {
    if (hash_comp(x, hashtab[hashval]) == 0) {
      hashtab[hashval] = NULL;
      Number_of_hash_entries--;
      return;
    }
    hash_value_rejections++;
    if (hash_value_rejections >= hashfillmax) {
      return;
    }
    hashval = ((hashval + 1) & (HASHSIZE));
  }
}


/* ----------------------------------------------------------------- */

/*
 * --------------------------------------------------------------------
 * mix -- mix 3 32-bit values reversibly.
 * For every delta with one or two bits set, and the deltas of all three
 * high bits or all three low bits, whether the original value of a,b,c
 * is almost all zero or is uniformly distributed,
 * If mix() is run forward or backward, at least 32 bits in a,b,c
 * have at least 1/4 probability of changing.
 * If mix() is run forward, every bit of c will change between 1/3 and
 * 2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
 * mix() takes 36 machine instructions, but only 18 cycles on a superscalar
 * machine (like a Pentium or a Sparc).  No faster mixer seems to work,
 * that's the result of my brute-force search.  There were about 2^^68
 * hashes to choose from.  I only tested about a billion of those.
 * --------------------------------------------------------------------
 */
#define mix(a, b, c) \
  { \
    a -= b; a -= c; a ^= (c >> 13); \
    b -= c; b -= a; b ^= (a << 8); \
    c -= a; c -= b; c ^= (b >> 13); \
    a -= b; a -= c; a ^= (c >> 12);  \
    b -= c; b -= a; b ^= (a << 16); \
    c -= a; c -= b; c ^= (b >> 5); \
    a -= b; a -= c; a ^= (c >> 3);  \
    b -= c; b -= a; b ^= (a << 10); \
    c -= a; c -= b; c ^= (b >> 15); \
  }

/*
 * --------------------------------------------------------------------
 * hash() -- hash a variable-length key into a 32-bit value
 * k       : the key (the unaligned variable-length array of bytes)
 * len     : the length of the key, counting by bytes
 * initval : can be any 4-byte value
 * Returns a 32-bit value.  Every bit of the key affects every bit of
 * the return value.  Every 1-bit and 2-bit delta achieves avalanche.
 * About 6*len+35 instructions.
 *
 * The best hash table sizes are powers of 2.  There is no need to do
 * mod a prime (mod is sooo slow!).  If you need less than 32 bits,
 * use a bitmask.  For example, if you need only 10 bits, do
 * h = (h & hashmask(10));
 * In which case, the hash table should have hashsize(10) elements.
 *
 * If you are hashing n strings (char **)k, do it like this:
 * for (i=0, h=0; i<n; ++i) h = hash( k[i], len[i], h);
 *
 * By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
 * code any way you wish, private, educational, or commercial.  It's free.
 *
 * See http://burtleburtle.net/bob/hash/evahash.html
 * Use for hash table lookup, or anything where one collision in 2^^32 is
 * acceptable.  Do NOT use for cryptographic purposes.
 * --------------------------------------------------------------------
 */
inline
PRIVATE unsigned
hash_f(void *x)
{
  register unsigned char  *k;           /* the key */
  register unsigned       length;       /* the length of the key */
  register unsigned       initval = 0;  /* the previous hash, or an arbitrary value */
  register unsigned       a, b, c, len;

  /* Set up the internal state */
  k   = ((hash_entry *)x)->structure;
  len = length = (unsigned)strlen(k);
  a   = b = 0x9e3779b9; /* the golden ratio; an arbitrary value */
  c   = initval;        /* the previous hash value */

  /*---------------------------------------- handle most of the key */
  while (len >= 12) {
    a += (k[0] + ((unsigned)k[1] << 8) + ((unsigned)k[2] << 16) + ((unsigned)k[3] << 24));
    b += (k[4] + ((unsigned)k[5] << 8) + ((unsigned)k[6] << 16) + ((unsigned)k[7] << 24));
    c += (k[8] + ((unsigned)k[9] << 8) + ((unsigned)k[10] << 16) + ((unsigned)k[11] << 24));
    mix(a, b, c);
    k   += 12;
    len -= 12;
  }

  /*------------------------------------- handle the last 11 bytes */
  c += length;
  switch (len) {
    /* all the case statements fall through */
    case 11:
      c += ((unsigned)k[10] << 24);
    case 10:
      c += ((unsigned)k[9] << 16);
    case 9:
      c += ((unsigned)k[8] << 8);
    /* the first byte of c is reserved for the length */
    case 8:
      b += ((unsigned)k[7] << 24);
    case 7:
      b += ((unsigned)k[6] << 16);
    case 6:
      b += ((unsigned)k[5] << 8);
    case 5:
      b += k[4];
    case 4:
      a += ((unsigned)k[3] << 24);
    case 3:
      a += ((unsigned)k[2] << 16);
    case 2:
      a += ((unsigned)k[1] << 8);
    case 1:
      a += k[0];
      /* case 0: nothing left to add */
  }
  mix(a, b, c);
  /*-------------------------------------------- report the result */
  return c & HASHSIZE;
}
