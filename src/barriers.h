#ifndef BARRIERS_H
#define BARRIERS_H

#include <barrier_types.h>

unsigned long *make_truemin(loc_min *Lmin);


loc_min *barriers(barrier_options opt);


int print_results(loc_min         *LM,
                  unsigned long   *tm,
                  barrier_options *opt);


void ps_tree(loc_min        *LM,
             unsigned long  *tm);


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


void free_rates(unsigned long length_rates);


void print_rates(unsigned long        n,
                 barrier_options      *opt,
                 barriers_rates_type  rate_files);


map_struc get_mapstruc(char           *p,
                       loc_min        *LM,
                       unsigned long  *tm);


unsigned long *compute_connected_component_states(loc_min       *lmin,
                                                  unsigned long *truemin);


int
print_rna_barriers_output(loc_min         *Lmin,
                          unsigned long   *truemin,
                          barrier_options *opt,
                          unsigned long   *mfe_component_true_min_indices);


void print_rates_of_mfe_component(unsigned long       *mfe_component_true_min_indices,
                                  barrier_options     *opt,
                                  barriers_rates_type rate_files);


void ps_tree_mfe_component(loc_min        *LM,
                           unsigned long  *tm,
                           unsigned long  *mfe_component_true_min_indices);


#endif
