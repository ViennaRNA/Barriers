/* Last changed Time-stamp: <2017-11-23 17:27:02 mtw> */
/* main.c */

#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "barrier_types.h"
#include "utils.h"
#include "barriers.h"
#include "hash_util.h"
#include "cmdline.h"

/* PRIVATE FUNCTIONS */
static barrier_options            opt;
static char                       *GRAPH;

static struct gengetopt_args_info args_info;
static int decode_switches(int  argc,
                           char **argv);


static void cleanup(char *,
                    loc_min *,
                    unsigned long *);


extern int  cut_point;
extern int  MYTURN;
extern char *tokenize(char *line);


extern char *costring(char *line);


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
  int           c;
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
      sec_structure     = (char *)calloc(len + 1, sizeof(char));
      protein_sequence  = (char *)calloc(len + 1, sizeof(char));
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

  LM = barriers(opt);
  if (opt.INFILE != stdin)
    fclose(opt.INFILE);

  tm = make_truemin(LM);

  if (opt.poset)
    mark_global(LM);

  if (cut_point > -1)
    opt.seq = costring(opt.seq);

  unsigned long *mfe_component_true_min_indices = NULL;
  if (opt.want_connected) {
    mfe_component_true_min_indices = compute_connected_component_states(LM, tm);
    fprintf(stderr, "computed mfe component\n");
    print_rna_barriers_output(LM, tm, &opt, mfe_component_true_min_indices);
  } else {
    c = print_results(LM, tm, &opt);
  }

  fflush(stdout);

  if (!opt.want_quiet) {
    if (opt.want_connected)
      ps_tree_mfe_component(LM, tm, 0, mfe_component_true_min_indices);
    else
      ps_tree(LM, tm, 0);
  }

  fprintf(stderr, "want_connected is %d\n", opt.want_connected);
  if (opt.want_connected && c == 0) {
    fprintf(stderr, "WARNING: landscape is not connected, skipping rates computation\n");
    cleanup(opt.seq, LM, tm);
    exit(102);
  }

  if (opt.rates != Barriers_no_rates || opt.microrates) {
    compute_rates(tm, opt.seq);
    if (!opt.want_quiet) {
      if (opt.want_connected)
        ps_tree_mfe_component(LM, tm, 1, mfe_component_true_min_indices);
      else
        ps_tree(LM, tm, 1);
    }

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
    unsigned long L1, L2;
    sscanf(args_info.path_arg[i], "%ld=%ld", &L1, &L2);
    if ((L1 > 0) && (L2 > 0)) {
      FILE        *PATH = NULL;
      char        tmp[30];
      path_entry  *path;

      path = backtrack_path(L1, L2, LM, tm);
      (void)sprintf(tmp, "path.%03ld.%03ld.txt", L1, L2);

      PATH = fopen(tmp, "w");
      if (PATH == NULL)
        nrerror("couldn't open path file");

      print_path(PATH, path, tm);
      /* fprintf(stderr, "%llu %llu\n", 0, MAXIMUM);   */
      fclose(PATH);
      fprintf(stderr, "wrote file %s\n", tmp);
      free(path);
    }
  }

  if (args_info.mapstruc_given) {
    FILE  *MAPFIN = NULL, *MAPFOUT = NULL;
    char  *line   = NULL, *token = NULL, *fname = "mapstruc.out";

    MAPFIN = fopen(args_info.mapstruc_arg, "r");
    if (MAPFIN == NULL)
      nrerror("couldn't open mapfile for reading");

    MAPFOUT = fopen(fname, "w");
    if (!MAPFOUT) {
      fprintf(stderr, "could not open mapstruc file %s for output\n", fname);
      errorcode = 101;
    }

    while (line = get_line(MAPFIN)) {
      map_struc myms;
      token = strtok(line, " \t");
      myms  = get_mapstruc(token, LM, tm);
      fprintf(MAPFOUT, "%s %6ld %6.2f %3ld %3ld %3ld %3ld\n", myms.structure, myms.n, myms.energy,
              myms.min, myms.truemin, myms.gradmin, myms.truegradmin);
      free(myms.structure);
      free(line);
    }
    fclose(MAPFIN);
    fclose(MAPFOUT);
  }

  free(mfe_component_true_min_indices);
  cleanup(opt.seq, LM, tm);
  exit(errorcode);
}


static void
cleanup(char          *seq,
        loc_min       *L,
        unsigned long *t)
{
  /* memory cleanup */
  free(seq);
  free(L);
  free(t);
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

  opt.max_print       = args_info.max_arg;
  opt.minh            = args_info.minh_arg;
  opt.poset           = args_info.poset_arg;
  opt.want_quiet      = args_info.quiet_given;
  opt.want_verbose    = args_info.verbose_given;
  opt.want_connected  = args_info.connected_given;
  opt.bsize           = args_info.bsize_given;
  opt.ssize           = args_info.ssize_given;
  opt.print_saddles   = args_info.saddle_given;
  opt.rates           = Barriers_no_rates;
  if (args_info.rates_given)
    opt.rates |= Barriers_both_rates; // print both files

  if (args_info.rates_text_file_given)
    opt.rates |= Barriers_text_rates; // print only text file

  opt.text_rates_file = args_info.rates_text_file_arg;
  if (args_info.rates_binary_file_given)
    opt.rates |= Barriers_binary_rates; // print only binary

  opt.binary_rates_file = args_info.rates_binary_file_arg;

  opt.microrates  = args_info.microrates_given;
  GRAPH           = args_info.graph_arg;
  opt.noLP_rate   = (args_info.noLP_rate_given) ? args_info.noLP_rate_arg : 1.;
  if (args_info.moves_given)
    opt.MOVESET = args_info.moves_arg;

  if (args_info.temp_given)
    opt.kT = args_info.temp_arg;

  for (i = 0; i < args_info.path_given; ++i) {
    unsigned long L1, L2;
    if (sscanf(args_info.path_arg[i], "%ld=%ld", &L1, &L2) != 2)
      nrerror("specifiy paths as e.g.  -P 1=3");
  }
  if (args_info.inputs_num > 1)
    nrerror("only one input file allowed");

  return 0;
}
