/* Last changed Time-stamp: <2001-03-08 17:37:43 ivo> */
/* hash_util.c */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "utils.h"
#include "hash_util.h"

#define   PRIVATE   static
#define   PUBLIC

static char UNUSED rcsid[] = "$Id: hash_util.c,v 1.1 2001/04/05 08:00:57 ivo Exp $";

/* modify hash_f(), hash_comp() and the typedef of hash_entry in hash_utils.h
   to suit your application */

PUBLIC void * lookup_hash (void *x);
PUBLIC int write_hash (void *x);
PUBLIC void delete_hash (void *x);
PUBLIC void kill_hash();
PUBLIC void initialize_hash();
PUBLIC int hash_comp(void *x, void *y);

PRIVATE unsigned hash_f (void *x);

#define HASHSIZE 33554432 -1  /* 2^25 -1   must be power of 2 -1 */
/* #define HASHSIZE 67108864 -1 */ /* 2^26 -1   must be power of 2 -1 */
/* #define HASHSIZE 16777216 -1 */ /* 2^24 -1   must be power of 2 -1 */ 
/* #define HASHSIZE 4194304 -1  */ /* 2^22 -1   must be power of 2 -1 */

PRIVATE void *hashtab[HASHSIZE+1];

PUBLIC unsigned long collisions=0;

/* ----------------------------------------------------------------- */

/* stolen from perl source */
char coeff[] = {
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1,
                61,59,53,47,43,41,37,31,29,23,17,13,11,7,3,1};

/* key must not be longer than 128 */
#pragma inline (hash_f)
PUBLIC unsigned hash_f(void *x)
{ 
  register char *s;
  register int i;
  register unsigned hash;

  s = ((hash_entry *)x)->structure;

  for (i=0,    hash = 0;
      /* while */ *s;
         s++,           i++ , hash *= 5 ) {
        hash += *s * coeff[i];
    }
  
  return ((hash) & (HASHSIZE)); /*divide through HASHSIZE for normalization */
}

PUBLIC int hash_comp(void *x, void *y) {
  return strcmp(((hash_entry *)x)->structure, ((hash_entry *)y)->structure);
}

/* ----------------------------------------------------------------- */
 
PUBLIC void * lookup_hash (void *x)  /* returns NULL unless x is in the hash */ 
{ 
  int hashval;

  hashval=hash_f(x);
  if (hashtab[hashval]==NULL) return NULL; 
  while (hashtab[hashval]){
    if (hash_comp(x,hashtab[hashval])==0) return hashtab[hashval];
    hashval = ((++hashval) & (HASHSIZE));
  }
  return NULL;
}

/* ----------------------------------------------------------------- */
    
PUBLIC int write_hash (void *x)   /* returns 1 if x already was in the hash */ 
{
  int hashval;
  
  hashval=hash_f(x);
  while (hashtab[hashval]){
    if (hash_comp(x,hashtab[hashval])==0) return 1;
    hashval = ((++hashval) & (HASHSIZE));
    collisions++;
  }
  hashtab[hashval]=x;
  return 0;
}
/* ----------------------------------------------------------------- */
  
PUBLIC void initialize_hash ()
{

}

/* ----------------------------------------------------------------- */

PUBLIC void kill_hash ()
{
  int i;
  
  for (i=0;i<HASHSIZE+1;i++) {
    if (hashtab[i]) {
      free (((hash_entry *)hashtab[i])->structure);
      free (hashtab[i]);
    }
    hashtab[i]=NULL;
  }
}

/* ----------------------------------------------------------------- */

PUBLIC void delete_hash (void *x)  /* doesn't work in case of collsions */
{                                  /* doesn't free anything ! */
  int hashval;
  
  hashval=hash_f(x);
  while (hashtab[hashval]){
    if (hash_comp(x,hashtab[hashval])==0) {
      hashtab[hashval]=NULL;
      return;
    }
    hashval = ((++hashval) & (HASHSIZE));
  }
}
/* ----------------------------------------------------------------- */

