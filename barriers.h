extern int      *make_truemin(loc_min *Lmin);
extern loc_min  *barriers(barrier_options opt);
extern void print_results(loc_min *LM, int *tm);
extern void ps_tree(loc_min *LM, int *tm);
extern path_entry *backtrack_path(int l1, int l2, loc_min *LM, int *truemin);
