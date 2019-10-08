/* Last changed Time-stamp: <2002-09-06 12:30:49 ivo> */
/* hash_util.h */

#ifndef _hash_util_h
#define _hash_util_h

extern void *lookup_hash(void *x);


extern int write_hash(void *x);


extern void delete_hash(void *x);


extern void kill_hash();


extern void initialize_hash();


typedef struct _hash_entry {
  char                *structure;     /* my structure */
  float               energy;         /* my energy */
  int                 basin;          /* which basin do I belong to */
  int                 GradientBasin;  /* for Gradient Basins */
  int                 ccomp;          /* in which connected component am I */
  int                 n;              /* my index in energy sorted list */
  struct _hash_entry  *down;          /* pointer to lowest neighbor */
  int                 *POV;           /* for Posets only */
} hash_entry;

#endif

/* End of file */
