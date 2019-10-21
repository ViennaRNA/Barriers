/* simple_set.h */
/* Last changed Time-stamp: <2001-07-04 20:31:06 ivo> */
#ifndef SIMPLE_SET_H
#define SIMPLE_SET_H

typedef struct {
  unsigned long basin;
  hash_entry    *hp;
} basinT;

typedef struct set {
  unsigned long num_elem;
  unsigned long max_elem;
  size_t        elem_size;
  basinT        *data;
} Set;

extern Set *new_set(int elems);


extern int set_add(Set    *set,
                   basinT *data);


extern void set_kill(Set *set);


extern int set_merge(Set        *s1,
                     const Set  *s2);


extern int set_find(Set     *set,
                    basinT  *data);


#endif
/* End of file */
