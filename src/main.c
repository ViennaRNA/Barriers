/*
 * main.c
 */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "barrier_types.h"
#include "utils.h"
#include "barriers.h"
#include "ringlist.h"
#include "cmdline.h"
#include "hash_table_linear_probing_lists/hash_tables.h"
#if HAVE_SECIS_EXTENSION
#include "secis_neighbors.h"
#endif

/* PRIVATE FUNCTIONS */
static barrier_options            opt;
static char                       *GRAPH;

static struct gengetopt_args_info args_info;
static int
decode_switches(int   argc,
                char  **argv);


static void
cleanup(char *,
        loc_min *,
        unsigned long *,
        vrna_hash_table_t *hash_table);


static char *program_name;

/*============================*/
int
main(int  argc,
     char *argv[])
{
  int           tmp;
  char          *line;
  loc_min       *LM;
  unsigned long *tm;
  unsigned long i, errorcode = 0;
  char          signal[100] = "", what[100] = "", stuff[100] = "";

  /* Parse command line */
  program_name = argv[0];

  opt.kT      = -300;
  opt.MOVESET = "";
  opt.minh    = 0.0000001;
  opt.label   = 0; /* normally, use numbers for minima */
  GRAPH       = NULL;

  /* Try to parse head to determine graph-type */
  decode_switches(argc, argv);

  if (args_info.inputs_num > 0) {
    opt.INFILE = fopen(args_info.inputs[0], "r");
    if (opt.INFILE == NULL)
      nrerror("can't open file");
  } else {
    opt.INFILE = stdin;
  }

  line = get_line(opt.INFILE);
  if (line == NULL) {
    fprintf(stderr, "Error in input file\n");
    exit(123);
  }

  opt.seq = (char *)space(strlen(line) + 1);
  sscanf(line, "%s %d %99s %99s %99s", opt.seq, &tmp, signal, what, stuff);
  opt.seq = tokenize(opt.seq);
  if (cut_point > -1)
    MYTURN = 1;

  if (strcmp(stuff, "\0") != 0 && strncmp(what, "Q", 1) == 0) {
    /* lattice proteins*/
    memset(opt.seq, 0, strlen(line) + 1);
    strcpy(opt.seq, stuff);
  }

  if ((!opt.poset) && (strcmp(signal, "::") != 0)) {
    int r, dim;
    /* in this case we have a poset file !!!! */
    r = sscanf(signal, "P:%d", &dim);
    if (r < 1) {
      fprintf(stderr,
              "Warning: obscure headline in input file\n");
      dim = 0;
    }

    if (dim > 0)
      opt.poset = dim;
  }

  if (opt.poset) {
    /* in this case we have a poset file !!!! */
    fprintf(stderr,
            "!!! Input data are a poset with %d objective functions\n",
            opt.poset);
    /* we have a SECIS design file */
    if (((GRAPH != NULL) && (strstr(GRAPH, "SECIS") != NULL))
        || (strncmp(what, "SECIS", 5) == 0)) {
#if HAVE_SECIS_EXTENSION
      int   len, max_m, min_as;
      char  *sec_structure, *protein_sequence;

      if (sscanf(what, "SECIS,%d,%d", &max_m, &min_as) < 2) {
        fprintf(stderr,
                "Error in input format for SECIS design !"
                "expected format: SECIS,INT,INT\n"
                "got: `%s'",
                what);
        exit(EXIT_FAILURE);
      }

      free(line);
      line              = get_line(opt.INFILE);
      len               = strlen(line);
      sec_structure     = (char *)space((len + 1)* sizeof(char));
      protein_sequence  = (char *)space((len + 1)* sizeof(char));
      sscanf(line, "%s %s", sec_structure, protein_sequence);

      if (opt.want_verbose)
        fprintf(stderr,
                "\nGraph is SECIS design with the following parameters:\n"
                "Structure:   %s\n"
                "Constraints: %s\n"
                "Protein sequence: %s\n"
                "Max. number of mutations : %d\n"
                "Min. alignment score (aa): %d\n\n",
                sec_structure,
                opt.seq,
                protein_sequence,
                max_m,
                min_as);

      initialize_SECIS(opt.seq, sec_structure, protein_sequence,
                       max_m, min_as);

      free(sec_structure);
      free(protein_sequence);
#else
      fprintf(stderr,
              "You need to reconfigure barriers with the --with-secis"
              " option\nto use barriers SECIS design extension\n");
      exit(EXIT_FAILURE);
#endif
    }
  }

  free(line);

  if (GRAPH == NULL)
    if (strlen(what))
      GRAPH = what;

  if (GRAPH == NULL)
    GRAPH = "RNA";

  opt.GRAPH = GRAPH;

  vrna_callback_ht_compare_entries *compare_function = hash_comp;
  vrna_callback_ht_hash_function   *hash_function = hash_function_uint64;
  vrna_callback_ht_free_entry      *free_hash_entry = barriers_free_hash_entry;
  vrna_hash_table_t hash_table = vrna_ht_init(HASHBITS, compare_function, hash_function, free_hash_entry);

  LM = barriers(opt, &hash_table);
  if (opt.INFILE != stdin)
    fclose(opt.INFILE);

  tm = make_truemin(LM);

  if (opt.poset)
    mark_global(LM);

  if (cut_point > -1)
    opt.seq = costring(opt.seq);

  unsigned long *mfe_component_true_min_indices = NULL;
  unsigned long max_mfe_comp_index              = 0;

  if (opt.want_connected) {
    mfe_component_true_min_indices = compute_connected_component_states(LM, tm);
    print_rna_barriers_output(LM, tm, &opt, mfe_component_true_min_indices);
    for (max_mfe_comp_index = 0;
         mfe_component_true_min_indices[max_mfe_comp_index] != 0;
         max_mfe_comp_index++);
  } else {
      print_results(LM, tm, &opt);
  }

  fflush(stdout);

  if (!opt.want_quiet) {
    if (opt.want_connected)
      ps_tree_mfe_component(LM, tm, mfe_component_true_min_indices);
    else
      ps_tree(LM, tm);
  }

  if (opt.rates != Barriers_no_rates || opt.microrates) {
    compute_rates(tm, opt.seq, &hash_table);

    if (opt.want_connected) {
      print_rates_of_mfe_component(mfe_component_true_min_indices, &opt, opt.rates);
      free_rates(tm[0]);
    } else {
      print_rates(tm[0], &opt, opt.rates);
    }
  }

  if (opt.poset)
    mark_global(LM);

  for (i = 0; i < args_info.path_given; ++i) {
    unsigned long l1, l2, l1_index_cc, l2_index_cc;
    sscanf(args_info.path_arg[i], "%ld=%ld", &l1, &l2);
    if ((l1 > 0) && (l2 > 0)) {
      FILE        *PATH = NULL;
      char        tmp[30];
      path_entry  *path;

      if (opt.want_connected) {
        /* map minima indices to connected component output indices */
        if (l1 > max_mfe_comp_index || l2 > max_mfe_comp_index) {
          fprintf(stderr,
                  "Error: one of the path indices is not in the connected component! l1=%ld, l2=%ld, maximum=%ld\n",
                  l1,
                  l2,
                  max_mfe_comp_index);
          exit(EXIT_FAILURE);
        }

        l1_index_cc = mfe_component_true_min_indices[l1 - 1];
        l2_index_cc = mfe_component_true_min_indices[l2 - 1];
        path        = backtrack_path(l1_index_cc, l2_index_cc, LM, tm, &hash_table);
      } else {
        path = backtrack_path(l1, l2, LM, tm, &hash_table);
      }

      (void)sprintf(tmp, "path.%03ld.%03ld.txt", l1, l2);

      PATH = fopen(tmp, "w");
      if (PATH == NULL)
        nrerror("couldn't open path file");

      if (opt.want_connected)
        print_path(PATH, path, tm, mfe_component_true_min_indices);
      else
        print_path(PATH, path, tm, NULL);

      /* fprintf(stderr, "%llu %llu\n", 0, MAXIMUM);   */
      fclose(PATH);
      fprintf(stderr, "wrote file %s\n", tmp);
      free(path);
    }
  }

  if (args_info.mapstruc_given) {
    FILE  *MAPFIN = NULL, *MAPFOUT = NULL;
    char  *line = NULL, *token = NULL, *fname = args_info.mapstruc_output_arg;

    MAPFIN = fopen(args_info.mapstruc_arg, "r");
    if (MAPFIN == NULL)
      nrerror("couldn't open mapfile for reading");

    MAPFOUT = fopen(fname, "w");
    if (!MAPFOUT) {
      fprintf(stderr, "could not open mapstruc file %s for output\n", fname);
      errorcode = 101;
    }

    while ((line = get_line(MAPFIN))) {
      map_struc myms;
      token = strtok(line, " \t");
      myms  = get_mapstruc(token, LM, tm, &hash_table);
      if (myms.structure != NULL) {
        if (args_info.connected_flag) {
          unsigned long mfe_comp_index;

          for (mfe_comp_index = 0; mfe_component_true_min_indices[mfe_comp_index] != 0;
               mfe_comp_index++)
            if (myms.truegradmin == mfe_component_true_min_indices[mfe_comp_index])
              break;

          if (mfe_comp_index < max_mfe_comp_index) {
            fprintf(MAPFOUT,
                    "%s %6ld %6.2f %3ld %3ld %3ld %3ld\n",
                    myms.structure,
                    myms.n,
                    myms.energy,
                    myms.min,
                    myms.truemin,
                    myms.gradmin,
                    mfe_comp_index);
          } else {
            fprintf(MAPFOUT, "structure not in mfe component!\n");
          }
        } else {
          fprintf(MAPFOUT,
                  "%s %6ld %6.2f %3ld %3ld %3ld %3ld\n",
                  myms.structure,
                  myms.n,
                  myms.energy,
                  myms.min,
                  myms.truemin,
                  myms.gradmin,
                  myms.truegradmin);
        }

        free(myms.structure);
      } else {
        fprintf(MAPFOUT, "structure not in hash\n");
      }

      free(line);
    }
    fclose(MAPFIN);
    fclose(MAPFOUT);
  }

  free(mfe_component_true_min_indices);
  cleanup(opt.seq, LM, tm, &hash_table);
  exit(errorcode);
}


static void
cleanup(char          *seq,
        loc_min       *L,
        unsigned long *t,
        vrna_hash_table_t *hash_table)
{
  /* memory cleanup */
  free(seq);
  free(L);
  free(t);
  vrna_ht_clear(*hash_table);
#if WITH_DMALLOC
  kill_hash(); /* freeing the hash takes unacceptably long */
#endif
  cmdline_parser_free(&args_info);
}


static int
decode_switches(int   argc,
                char  **argv)
{
  unsigned long i;

  if (cmdline_parser(argc, argv, &args_info) != 0)
    exit(1);

  if (args_info.help_given)
    cmdline_parser_print_help();

  if (args_info.full_help_given)
    cmdline_parser_print_full_help();

  opt.max_print       = (unsigned long)args_info.max_arg;
  opt.minh            = args_info.minh_arg;
  opt.poset           = args_info.poset_arg;
  opt.want_quiet      = args_info.quiet_flag;
  opt.want_verbose    = args_info.verbose_flag;
  opt.want_connected  = args_info.connected_flag;
  opt.bsize           = args_info.bsize_flag;
  opt.ssize           = args_info.ssize_flag;
  opt.print_saddles   = args_info.saddle_flag;
  opt.rates           = Barriers_no_rates;
  if (args_info.rates_given)
    opt.rates |= Barriers_both_rates; // print both files

  if (args_info.rates_text_file_given)
    opt.rates |= Barriers_text_rates; // print only text file

  opt.text_rates_file = args_info.rates_text_file_arg;
  if (args_info.rates_binary_file_given)
    opt.rates |= Barriers_binary_rates; // print only binary

  opt.binary_rates_file = args_info.rates_binary_file_arg;

  opt.microrates  = args_info.microrates_flag;
  GRAPH           = args_info.graph_arg;
  opt.noLP_rate   = (args_info.noLP_rate_given) ? args_info.noLP_rate_arg : 1.;
  if (args_info.moves_given)
    opt.MOVESET = args_info.moves_arg;

  if (args_info.temp_given)
    opt.kT = args_info.temp_arg;

  for (i = 0; i < args_info.path_given; ++i) {
    unsigned long L1, L2;
    if (sscanf(args_info.path_arg[i], "%ld=%ld", &L1, &L2) != 2)
      nrerror("specify paths as e.g.  -P 1=3");
  }
  if (args_info.inputs_num > 1)
    nrerror("only one input file allowed");

  return 0;
}
