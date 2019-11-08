#ifndef BARRIERS_H
#define BARRIERS_H

#include <stdbool.h>
#include "barrier_types.h"
#include "hash_table_linear_probing_lists/hash_tables.h"

void
set_barrier_options(barrier_options opt);


loc_min *
barriers(barrier_options opt, vrna_hash_table_t* hash_table);


unsigned long *
make_truemin(loc_min *Lmin);


void
check_neighbors(vrna_hash_table_t* hash_table);


void
mark_global(loc_min *Lmin);


void
print_results(loc_min         *Lmin,
              unsigned long   *truemin,
              barrier_options *opt);


void
print_rna_barriers_output(loc_min         *Lmin,
                          unsigned long   *truemin,
                          barrier_options *opt,
                          unsigned long   *mfe_component_true_min_indices);


char *
strip(char *s);


bool
is_bound(char *s);


unsigned long *
compute_connected_component_states(loc_min        *lmin,
                                   unsigned long  *truemin);


void
ps_tree(loc_min       *Lmin,
        unsigned long *truemin);


char *
get_taxon_label(int i);


path_entry *
backtrack_path(unsigned long  l1,
               unsigned long  l2,
               loc_min        *LM,
               unsigned long  *truemin,
               vrna_hash_table_t *hash_table);


void
print_path(FILE           *PATH,
           path_entry     *path,
           unsigned long  *tm,
           unsigned long  *mfe_component_true_min_indices);


map_struc
get_mapstruc(char           *p,
             loc_min        *LM,
             unsigned long  *tm,
             vrna_hash_table_t* hash_table);


void
print_rates(unsigned long       n,
            barrier_options     *opt,
            barriers_rates_type rate_files);


void
compute_rates(unsigned long *truemin,
              char          *sequence,
              vrna_hash_table_t* hash_table);


void
free_rates(unsigned long length_rates);


void
print_rates_of_mfe_component(unsigned long        *mfe_component_true_min_indices,
                             barrier_options      *opt,
                             barriers_rates_type  rate_files);


void
ps_tree_mfe_component(loc_min       *Lmin,
                      unsigned long *truemin,
                      unsigned long *mfe_component_true_min_indices);


#endif
