#ifndef BARRIERS_H
#define BARRIERS_H

#include <barrier_types.h>

int      *make_truemin(loc_min *Lmin);
loc_min  *barriers(barrier_options opt);
int print_results(loc_min *LM, int *tm, barrier_options *opt);
void ps_tree(loc_min *LM, int *tm, int rates);
path_entry *backtrack_path(int l1, int l2, loc_min *LM, int *truemin);
void print_path(FILE *PATH, path_entry *path, int *tm);
void mark_global(loc_min *Lmin);
void compute_rates(int *truemin, char *farbe);
void print_rates(int n, char *fname);
map_struc get_mapstruc(char *p, loc_min *LM, int *tm);

#endif
