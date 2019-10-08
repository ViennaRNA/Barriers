#ifndef BARRIERS_H
#define BARRIERS_H

#include <barrier_types.h>

unsigned long *make_truemin(loc_min *Lmin);


loc_min *barriers(barrier_options opt);


int print_results(loc_min         *LM,
                  unsigned long   *tm,
                  barrier_options *opt);


void ps_tree(loc_min        *LM,
             unsigned long  *tm,
             int            rates);


path_entry *backtrack_path(unsigned long  l1,
                           unsigned long  l2,
                           loc_min        *LM,
                           unsigned long  *truemin);


void print_path(FILE          *PATH,
                path_entry    *path,
                unsigned long *tm);


void mark_global(loc_min *Lmin);


void compute_rates(unsigned long  *truemin,
                   char           *farbe);


void print_rates(unsigned long  n,
                 char           *fname);


map_struc get_mapstruc(char           *p,
                       loc_min        *LM,
                       unsigned long  *tm);


#endif
