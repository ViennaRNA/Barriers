/* simple_set.h */
/* Last changed Time-stamp: <2001-07-04 20:31:06 ivo> */
#ifndef BARRIERS_SIMPLE_SET_H
#define BARRIERS_SIMPLE_SET_H

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

Set *new_set(int elems);


int set_add(Set    *set,
            basinT *data);


void set_kill(Set *set);


int set_merge(Set        *s1,
              const Set  *s2);


int set_find(Set     *set,
             basinT  *data);


#endif
/* End of file */
