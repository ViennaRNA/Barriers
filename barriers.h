extern int      *make_truemin(loc_min *Lmin);
extern loc_min  *barriers(barrier_options opt);
extern void print_results(loc_min *LM, int *tm, char *farbe);
extern void ps_tree(loc_min *LM, int *tm, int rates);
extern path_entry *backtrack_path(int l1, int l2, loc_min *LM, int *truemin);
extern void print_path(FILE *PATH, path_entry *path, int *tm);
extern void mark_global(loc_min *Lmin);
extern void compute_rates(int *truemin, char *farbe);
extern void print_rates(int n, char *fname);
