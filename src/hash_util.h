/* Last changed Time-stamp: <2002-09-06 12:30:49 ivo> */
/* hash_util.h */

#ifndef BARRIERS_HASH_UTIL_H
#define BARRIERS_HASH_UTIL_H

void *lookup_hash(void *x);


int write_hash(void *x);


void delete_hash(void *x);


void kill_hash();


void initialize_hash();


typedef struct _hash_entry {
  char                *structure;     /* my structure */
  float               energy;         /* my energy */
  unsigned long       basin;          /* which basin do I belong to */
  unsigned long       GradientBasin;  /* for Gradient Basins */
  unsigned long       ccomp;          /* in which connected component am I */
  unsigned long       n;              /* my index in energy sorted list */
  struct _hash_entry  *down;          /* pointer to lowest neighbor */
  int                 *POV;           /* for Posets only */
} hash_entry;

#endif

/* End of file */
