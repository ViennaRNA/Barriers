/* Last changed Time-stamp: <2003-11-28 13:22:18 mtw> */
/* main.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
/* #include <sys/types.h> */
#include <math.h>
#include <getopt.h>
#include "config.h"
#include "barrier_types.h"
#include "utils.h"
#include "barriers.h"
#include "hash_util.h"
         
/* PRIVATE FUNCTIONS */
static char UNUSED rcsid[] = "$Id: main.c,v 1.17 2003/11/28 12:28:26 mtw Exp $";
static void usage(int status);
static barrier_options opt;
static  char *GRAPH;
static int L1 = -1;
static int L2 = -1;

static struct option const long_options[] =
{
  {"quiet", no_argument, 0, 'q'},
  {"silent", no_argument, 0, 'q'},
  {"verbose", no_argument, 0, 'v'},
  {"help", no_argument, 0, 'h'},
  {"version", no_argument, 0, 'V'},
  {"label", no_argument, &opt.label, 1},
  {"max", required_argument, 0, 0},
  {"minh", required_argument, 0, 0},
  {"bsize", no_argument, &opt.bsize, 1},
  {"ssize", no_argument, &opt.ssize, 1},
  {"saddle", no_argument, &opt.print_saddles, 1},
  {"rates", no_argument, &opt.rates, 1},
  {"poset", required_argument, 0, 0 },
  {NULL, 0, NULL, 0}
};

static int decode_switches (int argc, char **argv);
static char* program_name;
/*============================*/
int main (int argc, char *argv[]) {
  int tmp;
  char *line;
  loc_min *LM;
  int *tm;
  int i;
  char signal[100]="", what[100]="", stuff[100]="";

  /* Parse command line */
  program_name = argv[0];

  opt.kT = -300;
  opt.MOVESET = "";
  opt.minh = 0.0000001;
  opt.label = 0; /* normally, use numbers for minima */
  GRAPH   = NULL;
    
  /* Try to parse head to determine graph-type */  
  i = decode_switches (argc, argv);

  if (i==argc-1) {
    opt.INFILE = fopen(argv[i], "r");
    if (opt.INFILE==NULL) nrerror("can't open file");
  } else {
    if (i<argc) usage(EXIT_FAILURE);
    opt.INFILE = stdin;
  }
  

  line = get_line(opt.INFILE);
  if (line == NULL) {
    fprintf(stderr,"Error in input file\n");
    exit(123);
  }
  opt.seq = (char *) space(strlen(line) + 1);
  sscanf(line,"%s %d %99s %99s %99s", opt.seq, &tmp, signal, what, stuff);
  if(strcmp(stuff, "\0")!=0 && strncmp(what, "Q", 1)==0){ /* lattice proteins*/ 
    memset(opt.seq, 0, strlen(line)+1);
    strcpy(opt.seq, stuff);
  }

  if ((!opt.poset)&&(strcmp(signal,"::")!=0)) {
    int r, dim;
    /* in this case we have a poset file !!!! */
    r=sscanf(signal,"P:%d",&dim);
    if (r<1) {
      fprintf(stderr,
	      "Warning: obscure headline in input file\n");
      dim = 0;
    }
    if (dim>0) opt.poset  = dim; 
  }

  if (opt.poset)
    fprintf(stderr,
	    "!!! Input data are a poset with %d objective functions\n",
	    opt.poset);
    
  free(line);

  if (GRAPH==NULL)
    if (strlen(what)) GRAPH = what;

  if (GRAPH==NULL) GRAPH="RNA";
  opt.GRAPH=GRAPH;
  
  LM = barriers(opt);
  if (opt.INFILE != stdin) fclose(opt.INFILE);
  tm = make_truemin(LM);

  print_results(LM,tm,opt.seq);
  fflush(stdout);

  if (!opt.want_quiet) ps_tree(LM,tm,0);
  
  if (opt.rates) {
    compute_rates(tm);
    if (!opt.want_quiet) ps_tree(LM,tm,1);
    print_rates(tm[0], "rates.out");
  }
  if (opt.poset) mark_global(LM);

  if ((L1>0) && (L2>0)) {
    FILE *PATH = NULL;
    char tmp[30];
    path_entry *path;

    path = backtrack_path(L1, L2, LM, tm);
    (void) sprintf(tmp, "path.%03d.%03d.txt", L1, L2);

    PATH = fopen (tmp, "w");
    if (PATH == NULL) nrerror("couldn't open path file");
    print_path(PATH, path, tm);
    /* fprintf(stderr, "%llu %llu\n", 0, MAXIMUM);   */
    fclose (PATH);
    fprintf (stderr, "wrote file %s\n", tmp);
  }

  /* memory cleanup */
  free(opt.seq);
  free(LM);
  free(tm);
#if WITH_DMALLOC
  kill_hash(); /* freeing the hash takes unacceptably long */
#endif
  exit(0);
}

static int decode_switches (int argc, char **argv)
{
  int c;
  int option_index = 0;

  while ((c = getopt_long (argc, argv, 
                           "q"  /* quiet or silent */
                           "v"  /* verbose */
                           "h"  /* help */
                           "V"  /* version */
			   "G:" /* GRAPH */
			   "M:" /* Move set */
			   "P:" /* backtrack path */
			   "T:",/* temperature for partition funktions */
                           long_options, &option_index)) != EOF)
    {
      switch (c)
        {
	case 0:
	  if (strcmp(long_options[option_index].name,"max")==0)
	    if (sscanf(optarg, "%d", &opt.max_print) == 0)
	      usage(EXIT_FAILURE);
	  if (strcmp(long_options[option_index].name,"minh")==0)
	    if (sscanf(optarg, "%lf", &opt.minh) == 0)
	      usage(EXIT_FAILURE);
	  if (strcmp(long_options[option_index].name,"poset")==0)
	    if (sscanf(optarg, "%d", &opt.poset) == 0)
	      usage(EXIT_FAILURE);
	  break;
        case 'q':               /* --quiet, --silent */
          opt.want_quiet = 1;
          break;
        case 'v':               /* --verbose */
          opt.want_verbose = 1;
          break;
        case 'V':
          printf ("barriers %s\n", VERSION);
          exit (0);
	  
        case 'h':
          usage (0);

	case 'G' :
	  GRAPH = optarg;
	  break;

	case 'M' :
	  opt.MOVESET = optarg;
	  break;
	case 'T':
          if (sscanf(optarg, "%lf", &opt.kT) == 0) usage(EXIT_FAILURE);
        break;

	case 'P':
	  if (sscanf(optarg, "%d=%d", &L1, &L2) != 2) usage(EXIT_FAILURE);
	  break;
        default:
          usage (EXIT_FAILURE);
        }
    }

  return optind;
}


/*==============================*/
static void usage(int status) {

  printf("%s - Compute local minima and energy barriers of landscape\n",
	  program_name);

  printf("Usage: %s [OPTION]... [FILE]\n", program_name);
  printf(
	 "Options:\n"
	 "-q, --quiet, --silent      be quiet, inhibit PS output\n"
	 "--verbose                  print more information\n"
	 "-h, --help                 display this help and exit\n"
	 "-V, --version              output version information and exit\n"
	 "-G <Graph>        define graph type.\n"
	 "-M Move-Set       select move-set\n"
	 "--bsize           log the basin sizes\n"
	 "--ssize           print out the saddle component sizes\n"
	 "--max <digit>     compute only the lowest <digit> local minima\n"
	 "--minh <de>       print only minima with barrier > de\n" 
	 "--saddle          log the saddle point structures\n"
	 "--rates           compute rates between macro states (basins)\n"
	 "-P <l1>=<l2>      backtrack path between lmins l2 and l1 (l1 < l2)\n"
         "--poset <n>       input is a poset from n objective functions\n"
	 );
  printf("\nFILE  must have RNAsubopt output-format sorted by energy\n\n");
  printf("Graph Types (-G graph) and Move Sets (-M mset)are:\n"
	 "  RNA             RNA secondary structures\n"
	 "  RNA-noLP        canonical RNA structures\n"
	 "      [no]Shift       with/out shift moves [default with]\n"
	 "  Q2              Spin Glass\n"
	 "      p               point mutation (default)\n"
	 "      c               flip 2nd half\n" 
	 "  Qa,ALPHA        a-letter Hamming graph.  \n"
	 "                      Specification of the ALPHAbet is optional\n"
	 "                      Default is 'ABC...'\n"
	 "  T               Phylogenetic Trees\n"
	 "      NNI             NNI moves [no other options yet]\n"
	 "  P               Permutations\n"
	 "      T               Transpositions [default]\n"
	 "      C               Canonical Transpositions\n"
	 "      R               Reversals\n"
	 "  X               Exchange Moves on balances +/- strings\n"
	 "  ?               General graph; adjacency list in file\n"
	 );
  
  exit (status);
}
