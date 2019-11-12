/* implement set as sorted array */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include "hash_tables.h"
#include "simple_set.h"


static int
comp_basinT(const void  *a,
            const void  *b)
{
  unsigned long A, B, index_a, index_b;

  A = ((basinT *)a)->basin;
  B = ((basinT *)b)->basin;
  if (A != B){
    if (A > B)
      return 1;
    else
      return -1;
  }

  if (((basinT *)a)->hp == NULL)
    return -1;

  if (((basinT *)b)->hp == NULL)
    return 1;

  /* else use energy or index in file */
  index_a = ((basinT *)a)->hp->n;
  index_b = ((basinT *)b)->hp->n;
  if(index_a > index_b)
    return 1;
  if(index_a < index_b)
    return -1;
  else
    return 0;
}


Set *
new_set(int elems)
{
  Set *set;

  set           = (Set *)space(sizeof(Set));
  set->data     = (basinT *)space(sizeof(basinT) * (elems + 1));
  set->max_elem = elems;
  return set;
}


int
set_find(Set    *set,
         basinT *data)
{
  int c = 1, i = 0, i2, i1 = 0;

  i2 = set->num_elem - 1;
  while (i1 <= i2) {
    i = (i1 + i2) / 2;
    c = comp_basinT(set->data + i, data);
    if (c > 0)
      i2 = i - 1;

    if (c < 0)
      i1 = i + 1;

    if (c == 0)
      break;
  }
  if (c == 0)
    return i;
  else
    return -i1 - 1;
}


int
set_add(Set     *set,
        basinT  *data)
{
  int pos;

  if ((pos = set_find(set, data)) >= 0)
    return 0;

  /* else insert before -pos-1 */
  pos = -pos - 1;
  set->num_elem++;
  if (set->max_elem <= set->num_elem) {
    set->max_elem *= 2;
    set->data     = (basinT *)xrealloc(set->data, sizeof(basinT) * (set->max_elem + 1));
  }

  memmove(set->data + pos + 1, set->data + pos, (set->num_elem - pos) * sizeof(basinT));
  set->data[pos] = *data;
  return 1;
}


void
set_kill(Set *set)
{
  free(set->data);
  free(set);
}


int
set_merge(Set       *s1,
          const Set *s2)
{
  int     i, i1, i2, num;
  basinT  *ndata;

  num = s1->num_elem + s2->num_elem + 1;
  if (num < s1->max_elem)
    num = s1->max_elem;

  ndata = (basinT *)space((num + 1) * sizeof(basinT));
  for (i = i1 = i2 = 0; i1 < s1->num_elem || i2 < s2->num_elem; i++) {
    int c;
    if (i1 == s1->num_elem)
      c = +1;
    else if (i2 == s2->num_elem)
      c = -1;
    else
      c = comp_basinT(s1->data + i1, s2->data + i2);

    ndata[i] = (c < 0) ?  s1->data[i1] : s2->data[i2];
    if (c <= 0)
      i1++;

    if (c >= 0)
      i2++;
  }
  free(s1->data);
  s1->data      = ndata;
  s1->max_elem  = num;
  s1->num_elem  = i;
  return i;
}


/* End of file */
