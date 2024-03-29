/*
 * barriers.c
 */

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#include "ringlist.h"
#include "stapel.h"
#include "utils.h"
#include "hash_tables.h"
#include "barrier_types.h"
#include "compress.h"
#include "treeplot.h"
#include "simple_set.h"
#include "rnamoves.h"   /* FK */
#if HAVE_SECIS_EXTENSION
#include "SECIS/secis_neighbors.h"
#endif

#include "barriers.h"

/* Tons of static arrays in this one! */
static char           *structure;  /* array for configuration */
static loc_min        *lmin;  /* array for local minima */

static double         **rate; /* rate matrix between basins */
static double         *dr;    /* increments to rate matrix  */

static unsigned long  n_lmin;
static unsigned long  max_lmin;
static unsigned long  n_saddle;
static double         minh;
static double         energy;
/* energy of last read structure (for check_neighbors) */
static int            *POV; /* list of last read POSET values */
static int            POV_size;
static double         mfe;  /* used for scaling Z */

static void           (*move_it)(char *);
static void           (*free_move_it)(void) = NULL;
static char           *(*pack_my_structure)(const char *);
static char           *(*unpack_my_structure)(const char *);

static double         kT = -1;

/* global switches */  /* defaults changed */
static int            print_saddles  = 1;
static int            bsize          = 1;
static int            shut_up        = 0;
static int            verbose        = 0;
static unsigned long  max_print      = 0;
static int            want_connected = 0;
static int            IS_RNA         = 0;
static int            print_labels   = 0;
static int            IS_arbitrary   = 0;

static int            maxlabellength = 0;

struct comp {
  Set           *basins;  /* set of basins connected by these saddles */
  char          *saddle;  /* one representative (first found) */
  unsigned long size;
};

static unsigned long  *truecomp;
static struct comp    *comp;
static unsigned long  max_comp = 1024, n_comp;
static int            do_rates      = 0;
static int            do_microrates = 0;
static double         noLP_rate     = 1.;
static int            ligand        = 0;

#define HASHSIZE (((unsigned long)1 << HASHBITS) - 1)
static hash_entry     *hpool;

/* private functions */
static void
Sorry(char *GRAPH);


static int
read_data(barrier_options opt,
          double          *energy,
          char            *strucb,
          int             len,
          int             *POV);


static void
merge_basins();


static int
path_cmp(const void *a,
         const void *b);


static void
backtrack_path_rec(unsigned long  l1,
                   unsigned long  l2,
                   const char     *tag,
                   hash_table_t* hash_table);


static void
walk_limb(hash_entry    *hp,
          unsigned long LM,
          int           inc,
          const char    *tag,
          hash_table_t *hash_table);


static void
merge_components(unsigned long  c1,
                 unsigned long  c2);


static int
comp_comps(const void *A,
           const void *B);


/* ----------------------------------------------------------- */

void
set_barrier_options(barrier_options opt)
{
  int isRNA2 = 0;

  print_saddles  = opt.print_saddles;
  bsize          = opt.bsize;
  shut_up        = opt.want_quiet;
  max_print      = opt.max_print;
  want_connected = opt.want_connected;
  minh           = opt.minh;
  verbose        = opt.want_verbose;
  print_labels   = opt.label;
  switch (opt.GRAPH[0]) {
    case 'R':   /* RNA secondary Structures */
      if (strncmp(opt.GRAPH, "RNA", 3) == 0) {
        int nolp = 0, shift = 0, i = 0;
        IS_RNA = 1;
        if (opt.kT <= -300)
          opt.kT = 37;

        kT      = 0.00198717 * (273.15 + opt.kT);     /* kT at 37C in kcal/mol */
        isRNA2  = !strncmp(opt.GRAPH, "RNA2", 4);     /* FK */

        if (isRNA2) {
          move_it       = RNA2_move_it;
          free_move_it  = RNA2_free;
        } else {
          move_it       = RNA_move_it;
          free_move_it  = RNA_free_rl;
        }

        pack_my_structure   = pack_structure;
        unpack_my_structure = unpack_structure;
        if (strstr(opt.GRAPH, "noLP")) {
          nolp      = 1;
          noLP_rate = opt.noLP_rate;
        }

        if (strlen(opt.MOVESET) > 0) {
          switch (opt.MOVESET[0]) {
            case 'S': /* Shift moves */
              if (strncmp(opt.MOVESET, "Shift", 5) == 0)
                shift = 1;

              break;

            case 's': /* Shift moves */
              if (strncmp(opt.MOVESET, "shift", 5) == 0)
                shift = 1;

              break;

            case 'n': /* no Shift moves */
              if (strncmp(opt.MOVESET, "noShift", 7) == 0 ||
                  strncmp(opt.MOVESET, "noshift", 7) == 0)
                shift = 0;

              break;

            case 'l': /* ligand */
              if (strncmp(opt.MOVESET, "ligand", 6) == 0)
                ligand = 1;

              /* move_it = RNA_move_it; */

              break;

            default:
              fprintf(stderr, "Unknown moveset %s\n", opt.MOVESET);
              exit(EXIT_FAILURE);
          }
        }

        for (i = 0; i < (int)strlen(opt.seq); i++)
          if (opt.seq[i] == 'T')
            opt.seq[i] = 'U';

        if (isRNA2)
          RNA2_init(opt.seq, shift, nolp);
        else
          RNA_init(opt.seq, shift, nolp);

        if (!shut_up)
          fprintf(stderr, "Graph is RNA with noLP=%d, Shift=%d, ligand=%d\n", nolp, shift, ligand);
      } else {
        Sorry(opt.GRAPH);
      }

      break;
    case 'Q':   /* Haming graphs */
      if (strcmp(opt.GRAPH, "Q2") == 0) {
        /* binary +- alphabet */
        if (strcmp(opt.MOVESET, "c") == 0) {
          move_it = SPIN_complement_move_it;
          if (!shut_up)
            fprintf(stderr, "Graph is Q2 with complementation moves\n");
        } else {
          move_it = SPIN_move_it;
          if (!shut_up)
            fprintf(stderr, "Graph is Q2\n");
        }

        pack_my_structure   = pack_spin;
        unpack_my_structure = unpack_spin;
      } else {
        int   alphabetsize = 0;
        int   numconv, i;
        char  *ALPHA;
        ALPHA   = (char *)space(strlen(opt.GRAPH) * sizeof(char));
        numconv = sscanf(opt.GRAPH, "Q%d,%s", &alphabetsize, ALPHA);
        switch (numconv) {
          case 2:
            if ((int)strlen(ALPHA) != alphabetsize)
              Sorry(opt.GRAPH);

            break;
          case 1:
            if ((alphabetsize <= 0) || (alphabetsize > 26))
              Sorry(opt.GRAPH);

            free(ALPHA);
            ALPHA = (char *)space(sizeof(char) * (alphabetsize + 1));
            for (i = 0; i < alphabetsize; i++)
              ALPHA[i] = (char)65 + i;
            break;
          default:
            Sorry(opt.GRAPH);
            break;
        }
        String_set_alpha(ALPHA);
        if (strcmp(opt.MOVESET, "c") == 0) {
          initialize_crankshaft();
          move_it = String_move_it_crankshaft;
          if (!shut_up)
            fprintf(stderr, "Graph is Q%d with Alphabet '%s' with crankshaft moves\n",
                    alphabetsize, ALPHA);
        } else {
          move_it = String_move_it;
          if (!shut_up)
            fprintf(stderr, "Graph is Q%d with Alphabet '%s'\n",
                    alphabetsize, ALPHA);
        }

        if (alphabetsize < 7) {
          ini_pack_em(opt);
          pack_my_structure   = pack_em;
          unpack_my_structure = unpack_em;
        } else {
          pack_my_structure   = strdup;
          unpack_my_structure = strdup;
        }
      }

      break;
    case 'P':   /* Permutations */
      switch (*opt.MOVESET) {
        case 'R':
          move_it = Reversal_move_it;
          break;
        case 'C':
          move_it = CTranspos_move_it;
          break;
        case 'T':
        default:

          move_it = Transpos_move_it;
      }
      pack_my_structure   = strdup;
      unpack_my_structure = strdup;
      if (!shut_up)
        fprintf(stderr, "Graph is Permutations with moveset %c\n",
                *opt.MOVESET ? *opt.MOVESET : 'T');

      break;
    case 'S':   /* multi objective SECIS design */
#if HAVE_SECIS_EXTENSION
      move_it             = SECIS_move_it;
      pack_my_structure   = strdup;
      unpack_my_structure = strdup;
#else
      fprintf(stderr,
              "You need to reconfigure barriers with the --with-secis"
              " option\nto use barriers SECIS design extension\n");
      exit(EXIT_FAILURE);
#endif
      break;
    case 'T':   /* Phylogenetic Trees */
      move_it             = NNI_move_it;
      pack_my_structure   = strdup;
      unpack_my_structure = strdup;
      if (!shut_up)
        fprintf(stderr, "Graph is Trees with NNI moves\n");

      break;
    case 'X': /* Johnson graph J(n,n/2) = balanced +/- with exchange moves */
      move_it             = EXCH_move_it;
      pack_my_structure   = pack_spin;
      unpack_my_structure = unpack_spin;
      break;
    case '?': /* General graph; adjacency list on file */
      move_it             = LIST_move_it;
      pack_my_structure   = strdup;
      unpack_my_structure = strdup;
      IS_arbitrary        = 1;
      break;
    default:
      Sorry(opt.GRAPH);
  }
  if (kT < 0) {
    if (opt.kT <= -300)
      kT = 1;
    else
      kT = opt.kT;
  }

  do_rates = opt.rates;
  if (opt.microrates) {
    do_microrates = opt.microrates;
    do_rates      = opt.microrates;
  }
}


static void
Sorry(char *GRAPH)
{
  fprintf(stderr, "Graph \"%s\" is not implemented\n", GRAPH);
  exit(-2);
}


static FILE           *Mergefile  = NULL;
static unsigned long  Read_lines  = 0;
loc_min *
barriers(barrier_options opt, hash_table_t* hash_table)
{
  int     length;
  double  new_en = 0;

  fprintf(stderr, "hashbits = %u\n", HASHBITS);
  hpool = (hash_entry *)space((HASHSIZE + 1) * sizeof(hash_entry));
  set_barrier_options(opt);

  length    = (int)strlen(opt.seq);
  max_lmin  = 16383;
  lmin      = (loc_min *)space((max_lmin + 1) * sizeof(loc_min));
  n_lmin    = 0;

  structure      = (char *)space((length + 2) * sizeof(char));
  comp      = (struct comp *)space((max_comp + 1) * sizeof(struct comp));
  truecomp  = (unsigned long *)space((max_comp + 1) * sizeof(unsigned long));
  if (opt.poset) {
    POV_size  = opt.poset;
    POV       = (int *)space(sizeof(int) * opt.poset);
  } else {
    POV = NULL;
  }

  ini_stapel(length);
  if (opt.ssize) {
    Mergefile = fopen("saddles.txt", "w");
    if (!Mergefile)
      fprintf(stderr, "can't open saddle file\n");
  }

  while (read_data(opt, &new_en, structure, length, POV)) {
    if (Read_lines == 0)
      mfe = energy = new_en;

    if (new_en < energy)
      nrerror("unsorted list!\n");

    if (new_en > energy) {
      /* new energy band started */
      merge_basins();
      /* fprintf(stderr, "%d %d\n", readl, lmin[1].my_pool); */
      n_comp = 0;
    }

    energy = new_en;
    Read_lines++;
    move_it(structure);        /* generate all neighbor of configuration */
    /* fprintf(stderr, "M%s\n", form); */
    check_neighbors(hash_table);    /* flood the energy landscape */
    reset_stapel();
    if ((n_saddle + 1 == max_print) && (!opt.rates))
      break;  /* we've found all we want to know */
  }
  switch (opt.GRAPH[0]) {
    case 'Q':
      if (strcmp(opt.MOVESET, "c") == 0)
        Q_mem_cleanup();

      break;
    default:
      break;
  }
  merge_basins();
  if (Mergefile)
    fclose(Mergefile);

  if (!shut_up)
    fprintf(stderr,
            "read %ld structures, to find %ld saddles\n",
            Read_lines, n_saddle);

  if (max_print == 0 || max_print > n_lmin)
    max_print = n_lmin;

  lmin[0].fathers_pool  = n_lmin;  /* store size here; pfs 03 2001 */
  lmin[0].E_saddle      = energy + 0.001;
  lmin[0].energy        = lmin[1].energy;

  if (!do_rates) {
    if (free_move_it)
      free_move_it();

    free_stapel();
  }

  free(structure);
  fflush(stdout);
  if (!shut_up){
    //fprintf(stderr, "%lu hash table collisions\n", collisions);
    fprintf(stderr, "%lu hash table collisions\n", ht_collisions(*hash_table));
  }
  free(truecomp);
  free(comp);
  return lmin;
}


unsigned long *
make_truemin(loc_min *Lmin)
{
  unsigned long *truemin, nlmin, i, ii, f;
  int           *is_connected = NULL;     /* stores which mins are connected */

  nlmin   = Lmin[0].fathers_pool;
  truemin = (unsigned long *)space((nlmin + 1) * sizeof(unsigned long));
  /* truemin[0] = nlmin; */

  /* Ensure that max_print mins are printed when using the --connected option.
   * Initialize list to keep track of connected mins. */
  is_connected = (int*) space((n_lmin + 1) * sizeof(int));
  for (i = 0; i <= n_lmin; i++) { is_connected[i] =  0; } /* assume disconn */
  is_connected[1] = 1;                    /* min 1 is connected */

  for (ii = i = 1; (i <= max_print) && (ii <= n_lmin); ii++) {
    f = lmin[ii].father;
    if (!f)
      lmin[ii].E_saddle = energy + 0.000001;
    if (is_connected[f])  			/* i != 1 is connected iff its father is */
      is_connected[ii] = 1;

    /* A min is true iff 1) its barrier height is at least minh, and 2) if the
     * want_connected option is given, it has to be connected to the mfe. */
    if (((lmin[ii].E_saddle - lmin[ii].energy - (float)minh + FLT_EPSILON >= 0.)
          && (!want_connected || is_connected[ii]))
        || (nlmin == 1)) {
      truemin[ii] = i++;
    } else {
      /* ii is not a truemin */
      lmin[f].Z   += lmin[ii].Z;
      lmin[f].Zg  += lmin[ii].Zg;
    }
  }
  truemin[0] = i - 1;
  return truemin;
}


/*=============================================================*/

static int
read_data(barrier_options opt,
          double          *energy,
          char            *strucb,
          int             len,
          int             *POV)
{
  int     r, l;
  char    *line;
  char    *token;

#ifdef _DEBUG_POSET_
  static  count = 1;
#endif

  line = get_line(opt.INFILE);

  if (line == NULL)
    return 0;

  if (strlen(line) == 0)
    return 0;

  token = strtok(line, " \t");
  if (token == NULL)
    return 0;

  l = (int)strlen(token);
  if (l < 1)
    return 0;

  if (cut_point == -1) {
    strcpy(strucb, token);
  } else {
    char *mystr;
    mystr = (char *)space((l + 1) * sizeof(char));
    strcpy(mystr, token);
    mystr = tokenize(mystr);
    strcpy(strucb, mystr);
    free(mystr);
  }

  token = strtok(NULL, " \t");
  if (token == NULL) {
    fprintf(stderr, "Error in input file\n");
    exit(123);
  }

  r = sscanf(token, "%lf", energy);
  if (r < 1) {
    fprintf(stderr, "Error in input file\n");
    exit(124);
  }

  if (IS_arbitrary) {
    /* record the maximal length of token name for output formatting */
    if (l > maxlabellength)
      maxlabellength = l;

    if (l > len) {
      fprintf(stderr, "read_data():\n%s\n label too long !!\n", strucb);
      exit(111);
    }
  }

#if 0
  /*
   * removed because in the lattice protein case, the sequence is one
   * character longer than the actual SAWs
   */
  else {
    if (l != len) {
      fprintf(stderr, "read_data():\n%s\n unequal length !!\n", strucb);
      exit(112);
    }
  }
#endif

  if (opt.poset) {
    int i, x;
    for (i = 0; i < opt.poset; i++) {
      token = strtok(NULL, " \t");
      if (token == NULL) {
        fprintf(stderr, "Error in input file\n");
        exit(125);
      }

      r = sscanf(token, "%d", &x);
      if (r != 1) {
        fprintf(stderr, "Error in input file\n");
        exit(126);
      }

      POV[i] = x;
    }
#ifdef _DEBUG_POSET_
    {
      int i;
      fprintf(stderr, "POV[%4d] = {", count);
      for (i = 0; i < opt.poset; i++) {
        fprintf(stderr, "%2d", POV[i]);
        if (i < opt.poset - 1)
          fprintf(stderr, ",");
      }
      fprintf(stderr, "}\n");
    }
    count++;
#endif
  }

  if (IS_arbitrary) {
    token = strtok(NULL, " \t");
    if (token == NULL)
      put_ADJLIST(":");
    else
      put_ADJLIST(token);
  }

  /* the rest of the line is junk an can safely (?) be ignored */

  free(line);

  return 1;
}


/*======================*/
void
check_neighbors(hash_table_t* hash_table)
{
  char          *p, *pp, *pform;
  unsigned long basin, obasin = ULONG_MAX;
  hash_entry    *hp, h, *down = NULL;
  Set           *basins;
  basinT        b;

  double        Zi;
  unsigned long min_n   = ULONG_MAX;      /* index of lowest neighbor  */
  unsigned long gradmin = 0;              /* for Gradient Basins */
  unsigned long is_min  = 1;
  unsigned long ccomp   = 0;              /* which connected component */

  basins = new_set(10);

  Zi = exp((mfe - energy) / kT);

  /* foreach neighbor structure of configuration "Structure" */
  while ((p = pop())) {
    if (IS_RNA)
      if (p[strlen(p) - 1] == 'D')
        p[strlen(p) - 1] = '\0';

    pp          = pack_my_structure(p);
    h.structure = pp;


    //hp = lookup_hash(&h);
    hp = (hash_entry *)ht_get(*hash_table, &h);

    if (hp && POV_size) {
      /* need to check if h is dominated by hp */
      int i;
      for (i = 0; i < POV_size; i++)
        /* printf(" %d",hp->POV[i]); */
        if (POV[i] < hp->POV[i]) {
          hp = NULL;
          break;
        }
    }

    /* check whether we've seen the structure before */
    if (hp) {
      /*
       * because we've seen this structure before, it already
       * belongs to the basin of attraction of a local minimum
       */
      basin = hp->basin;
      if (hp->energy < energy)
        is_min = 0;                        /* should we use hp->n here? */

      if (hp->n < min_n) {
        /* find lowest energy neighbor */
        min_n   = hp->n;
        gradmin = hp->GradientBasin;
        down    = hp;
      }

      /* careful: if input has higher precision than FLT_EPSILON
       * bad things will happen */
      if (fabs(hp->energy - energy) <= FLT_EPSILON * fabs(energy)) {
        unsigned long tc;
        tc = hp->ccomp;
        while (tc != truecomp[tc])
          tc = truecomp[tc];
        if (ccomp == 0) {
          ccomp = tc;
        } else {
          while (ccomp != truecomp[ccomp])
            ccomp = truecomp[ccomp];
          if (ccomp != tc)
            merge_components(tc, ccomp);

          ccomp = truecomp[ccomp];
        }
      }

      /*
       * the basin of attraction of this local minimum may have been
       * merged with the basin of attraction of an energetically
       * "deeper" local minimum in a previous step
       * go and find this "deeper" local minimum!
       */
      while (lmin[basin].father)
        basin = lmin[basin].father;

      /* put the "deepest" local minimum into the basins-list */
      if (basin != obasin || obasin == ULONG_MAX) {
        b.hp    = hp;
        b.basin = basin;
        set_add(basins, &b);
      }

      obasin = basin;
    }

    free(pp);
  }

  /*
   * pack read structure from subopt for putting into the hash
   * fprintf(stderr,"F%s\n",form);
   */
  pform = pack_my_structure(structure);

  if (ccomp == 0) {
    /* new compnent */
    Set *set;
    set = new_set(10);
    if (++n_comp > max_comp) {
      max_comp  *= 2;
      comp      = (struct comp *)xrealloc(comp, (max_comp + 1) * sizeof(struct comp));
      truecomp  = (unsigned long *)xrealloc(truecomp, (max_comp + 1) * sizeof(unsigned long));
    }

    comp[n_comp].basins = set;
    comp[n_comp].saddle = pform;
    comp[n_comp].size   = 0;
    truecomp[n_comp]    = ccomp = n_comp;
  }

  if (is_min) {
    basinT b;
    /* Structure is a "new" local minimum */
    gradmin = ++n_lmin;        /* for Gradient Basins */
    down    = NULL;
    /* need to allocate more space for the lmin-list */
    if (n_lmin > max_lmin) {
      fprintf(stderr, "increasing lmin array to %ld\n", max_lmin * 2);
      lmin = (loc_min *)xrealloc(lmin, (max_lmin * 2 + 1) * sizeof(loc_min));
      memset(lmin + max_lmin + 1, 0, max_lmin);
      max_lmin *= 2;
    }

    /* store configuration "Structure" in lmin-list */
    lmin[n_lmin].father       = 0;
    lmin[n_lmin].structure    = pform;
    lmin[n_lmin].energy       = energy;
    lmin[n_lmin].my_GradPool  = 0;
    lmin[n_lmin].my_pool      = 1;
    lmin[n_lmin].Z            = Zi;
    lmin[n_lmin].Zg           = 0;
    b.basin                   = n_lmin;
    b.hp                      = NULL;
    set_add(basins, &b);
  } else {
    comp[ccomp].size++;
  }

  set_merge(comp[ccomp].basins, basins);

  {
    unsigned long i_lmin;
    i_lmin = (is_min) ? n_lmin : basins->data[0].basin;
    set_kill(basins);
    /* store configuration "Structure" in hash table */
    if (Read_lines > HASHSIZE) {
      fprintf(stderr,
              "Error: Structure in line %ld could not be written to the hash table! Please restrict the input or recompile with --with-hash-bits and a higher value.\n",
              Read_lines);
      exit(EXIT_FAILURE);
    }

    hp = hpool + Read_lines - 1;  /* (hash_entry *) space(sizeof(hash_entry)); */
    if (POV_size) {
      int i;
      hp->POV = (int *)space(sizeof(int) * POV_size);
      for (i = 0; i < POV_size; i++)
        hp->POV[i] = POV[i];
    }

    hp->structure     = pform;
    hp->energy        = energy;
    hp->basin         = i_lmin;
    hp->GradientBasin = gradmin;    /* for Gradient Basins */
    hp->down          = down;
    hp->ccomp         = ccomp;
    hp->n             = Read_lines;
    lmin[gradmin].my_GradPool++;
    lmin[gradmin].Zg += Zi;

    //int write_result = write_hash(hp);
    int write_result = ht_insert(*hash_table, hp);
    if (write_result == 1) {
      hash_entry *foo = NULL;
      //foo = (hash_entry *)lookup_hash(hp);
      foo = (hash_entry *)ht_get(*hash_table, hp);
      fprintf(stderr, "%s\n", unpack_my_structure(foo->structure));
      fprintf(stderr, "%s\n", unpack_my_structure(hp->structure));
      nrerror("duplicate structure");
    }

    if (write_result == -1) {
      if (!shut_up)
        fprintf(stderr, "%lu hash table collisions\n", ht_collisions(*hash_table));

      exit(EXIT_FAILURE);
    }
  }

  if ((is_min) && (POV_size))
    lmin[n_lmin].POV = hp->POV;
}


static void
merge_basins()
{
  unsigned long c, i, t;

  for (i = t = 1; i <= n_comp; i++) {
    if (truecomp[i] == i)
      comp[t++] = comp[i];
    else
      set_kill(comp[i].basins);
  }
  n_comp = t - 1;
  qsort(comp + 1, n_comp, sizeof(struct comp), comp_comps);
  for (c = 1; c <= n_comp; c++) {
    /*
     * foreach connected component
     * merge all lmins connected by this component
     */
    static unsigned long  false_lmin = 0;
    unsigned long         i, father, pool = 0;
    double                Z = 0;
    basinT                *basins;

    if (Mergefile && (comp[c].basins->num_elem > 1)) {
      const char  format[2][16] = {
        "%13.5f %4d %s", "%6.2f %4d %s"
      };
      char        *saddle;
      saddle = unpack_my_structure(comp[c].saddle);
      fprintf(Mergefile, format[IS_RNA], energy, comp[c].size, saddle);
      free(saddle);
      for (i = 0; i < comp[c].basins->num_elem; i++)
        fprintf(Mergefile, " %2ld", comp[c].basins->data[i].basin);
      fprintf(Mergefile, "\n");
    }

    basins  = comp[c].basins->data;
    father  = basins[0].basin;
    while (lmin[father].father)
      father = lmin[father].father;

    for (i = 1; i < comp[c].basins->num_elem; i++) {
      unsigned long ii, l, r;
      ii = basins[i].basin;
      while (lmin[ii].father)
        ii = lmin[ii].father;
      if (ii != father) {
        if (ii < father) {
          unsigned long tmp;
          tmp     = ii;
          ii      = father;
          father  = tmp;
          l       = 0;
          r       = i;
        } else {
          l = i;
          r = 0;
        }

        /* going to merge ii with father  */
        if ((!max_print) || (ii <= max_print + false_lmin)) {
          /* found the saddle for a basin we're gonna print */
          if (energy - lmin[ii].energy >= minh)
            n_saddle++;
          else
            false_lmin++;
        }

        lmin[ii].father   = father;
        lmin[ii].saddle   = comp[c].saddle;
        lmin[ii].E_saddle = energy;
        lmin[ii].left     = basins[l].hp;
        lmin[ii].right    = basins[r].hp;
        if (bsize) {
          lmin[ii].fathers_pool = lmin[father].my_pool;
          pool                  += lmin[ii].my_pool;
          Z                     += lmin[ii].Z;
        }
      }
    }
    if (bsize) {
      lmin[father].my_pool  += pool + comp[c].size;
      lmin[father].Z        += Z + comp[c].size * exp((mfe - energy) / kT);
    }

    set_kill(comp[c].basins);
  }
}


void
mark_global(loc_min *Lmin)
{
  unsigned long i, j, k;
  unsigned long n_gmin;
  loc_min       *G;

  G = (loc_min *)space(n_lmin * sizeof(loc_min));

  Lmin[1].global  = (char)1;
  G[1]            = Lmin[1];
  n_gmin          = 1;

  for (i = 1; i <= n_lmin; i++) {
    Lmin[i].global = (char)1;
    for (j = 1; j <= n_gmin; j++) {
      int dom;
      dom = 1;
      for (k = 0; k < POV_size; k++)
        if (G[j].POV[k] >= Lmin[i].POV[k])
          dom = 0;

      if (dom) {
        Lmin[i].global  = (char)0;
        j               = n_gmin + 1;
      }
    }
    if (Lmin[i].global) {
      n_gmin++;
      G[n_gmin] = Lmin[i];
    }
  }
}


/*====================*/
void
print_results(loc_min         *Lmin,
              unsigned long   *truemin,
              barrier_options *opt)
{
  char          *sequence = opt->seq;
  unsigned long i, ii, j, k, n;
  char          *struc = NULL, *laststruc = NULL;
  char          *format = NULL, *formatA = NULL, *formatB = NULL;
  bool          otherformat = false;

  if (POV_size)
    fprintf(stderr, " POV_size = %d\n", POV_size);

  if (IS_arbitrary) {
    char tfor[100];
    sprintf(tfor, "%%4d %%-%ds %%6.2f %%4d %%6.2f", maxlabellength);
    format = tfor;
  } else if (IS_RNA) {
    formatA = "%4d %s %6.2f %4d %6.2f";
    formatB = "%4d %s  %6.2f %4d %6.2f";
    format  = formatA;
  } else {
    format = "%4d %s %13.5f %4d %13.5f";
  }

  if (verbose)
    printf("Using output format string '%s'\n", format);

  n_lmin = Lmin[0].fathers_pool;

  printf("     %s\n", sequence);
  for (i = 1, k = 0; i <= n_lmin; i++, k++) {
    unsigned long f;
    if ((ii = truemin[i]) == 0)
      continue;

    struc = unpack_my_structure(Lmin[i].structure);
    if (cut_point > -1)
      struc = costring(struc);

    n = strlen(struc);
    f = Lmin[i].father;
    if (f > 0)
      f = truemin[f];

    if (POV_size) {
      int jj;
      printf("%4ld %s ", ii, struc);

      if (IS_RNA)
        printf("%6.2f ", Lmin[i].energy);
      else
        printf("%13.5f ", Lmin[i].energy);

      for (jj = 0; jj < POV_size; jj++)
        printf("%6d ", Lmin[i].POV[jj]);
      if (IS_RNA)
        printf("%4ld %6.2f", f,
               Lmin[i].E_saddle - Lmin[i].energy);
      else
        printf("%4ld %13.5f", f,
               Lmin[i].E_saddle - Lmin[i].energy);

      if (Lmin[i].global)
        printf(" *");
      else
        printf(" .");
    } else {
      if (IS_RNA && (ligand == 1)) {
        if (strstr(struc, "*") == NULL) {
          format      = formatB;
          otherformat = true;
        }
      }

      printf(format, ii, struc, Lmin[i].energy, f,
             Lmin[i].E_saddle - Lmin[i].energy);
      if (otherformat) {
        format      = formatA;
        otherformat = false;
      }
    }

    if (print_saddles) {
      if (Lmin[i].saddle) {
        char *saddlestruc = unpack_my_structure(Lmin[i].saddle);
        printf(" %s", saddlestruc);
        free(saddlestruc);
      } else {
        printf(" ");
        for (j = 0; j < n; j++)
          printf("~");
      }
    }

    if (bsize)
      printf(" %12ld %8ld %10.6f %8ld %10.6f",
             Lmin[i].my_pool, Lmin[i].fathers_pool, mfe - kT * log(lmin[i].Z),
             Lmin[i].my_GradPool, mfe - kT * log(lmin[i].Zg));

    printf("\n");

    if (laststruc != NULL)
      free(laststruc);

    laststruc = strdup(struc);
    free(struc);
  } /* end for */
  free(laststruc);
}


void
print_rna_barriers_output(loc_min         *Lmin,
                          unsigned long   *truemin,
                          barrier_options *opt,
                          unsigned long   *mfe_component_true_min_indices)
{
  char          *sequence = opt->seq;
  unsigned long i, ii, k, n;
  char          *struc = NULL, *laststruc = NULL;

  char          *format = NULL, *formatA = NULL, *formatB = NULL;
  bool          otherformat = false;

  if (IS_RNA) {
    formatA = "%4d %s %6.2f %4d %6.2f";
    formatB = "%4d %s  %6.2f %4d %6.2f";
    format  = formatA;
  }

  if (verbose)
    printf("Using output format string '%s'\n", format);

  n_lmin  = Lmin[0].fathers_pool;
  i       = 1;
  unsigned long max_mfe_comp_size = 0;
  if (mfe_component_true_min_indices != NULL) {
    while (mfe_component_true_min_indices[max_mfe_comp_size] != 0)
      max_mfe_comp_size++;
    //fprintf(stderr, "%ld states in mfe component!\n", max_mfe_comp_size);
  }

  printf("     %s\n", sequence);
  for (k = 0; i <= n_lmin; i++, k++) {
    unsigned long f;
    if ((ii = truemin[i]) == 0)
      continue;

    unsigned long j = 0;

    if (mfe_component_true_min_indices != NULL) {
      int is_in_mfe_comp = 0;
      for (j = 0; mfe_component_true_min_indices[j] != 0; j++) {
        if (mfe_component_true_min_indices[j] == ii) {
          is_in_mfe_comp = 1;
          break;
        }
      }
      if (!is_in_mfe_comp)
        continue;
    }

    struc = unpack_my_structure(Lmin[i].structure);
    if (cut_point > -1)
      struc = costring(struc);

    n = strlen(struc);
    f = Lmin[i].father;
    if (f > 0) {
      f = truemin[f];
      if (mfe_component_true_min_indices != NULL) {
        unsigned long mfe_comp_index = 0;
        while (mfe_component_true_min_indices[mfe_comp_index] != 0) {
          if (f == mfe_component_true_min_indices[mfe_comp_index])
            break;

          mfe_comp_index++;
        }
        if (mfe_comp_index >= max_mfe_comp_size)
          //fprintf(stderr, "Error: father %ld is not in mfe component!\n", f);
          continue;

        f = mfe_comp_index + 1;
      }
    }

    if (IS_RNA && (ligand == 1)) {
      if (strstr(struc, "*") == NULL) {
        format      = formatB;
        otherformat = true;
      }
    }

    // numbers with step size 1
    if (mfe_component_true_min_indices != NULL)
      ii = j + 1;

    printf(format, ii, struc, Lmin[i].energy, f,
           Lmin[i].E_saddle - Lmin[i].energy);
    if (otherformat) {
      format      = formatA;
      otherformat = false;
    }

    if (print_saddles) {
      if (Lmin[i].saddle) {
        char *saddlestruc = unpack_my_structure(Lmin[i].saddle);
        printf(" %s", saddlestruc);
        free(saddlestruc);
      } else {
        printf(" ");
        for (j = 0; j < n; j++)
          printf("~");
      }
    }

    if (bsize)
      printf(" %12ld %8ld %10.6f %8ld %10.6f",
             Lmin[i].my_pool, Lmin[i].fathers_pool, mfe - kT * log(lmin[i].Z),
             Lmin[i].my_GradPool, mfe - kT * log(lmin[i].Zg));

    printf("\n");

    if (laststruc != NULL)
      free(laststruc);

    laststruc = strdup(struc);
    free(struc);
  } /* end for */
  free(laststruc);
}


/*====================*/
/* remove additional characters from structure, such as
 * '*','A','B',... */
char *
strip(char *s)
{
  char *p = strdup(s);

  if (p[strlen(p) - 1] == '*') /* add more characters here is required */
    p[strlen(p) - 1] = '\0';

  return p;
}


/*====================*/
/* check if a structure is marked ligand/protein-bound, ie it has an
 * extra char (eg '*') appended */
bool
is_bound(char *s)
{
  bool val = false;

  if (s[strlen(s) - 1] == '*') /* add more characters here is required */
    val = true;

  return val;
}


unsigned long *
compute_connected_component_states(loc_min        *lmin,
                                   unsigned long  *truemin)
{
  unsigned long nlmin                 = lmin[0].fathers_pool;
  unsigned long *mfe_component_minima = (unsigned long *)space(sizeof(unsigned long) * (truemin[0] + 1));
  unsigned long ii;

  unsigned long  star_mfe_index    = ULONG_MAX;
  unsigned long  normal_mfe_index  = ULONG_MAX;
  float         star_mfe          = FLT_MAX;
  float         normal_mfe        = FLT_MAX;

  loc_min       *current_lm;
  char          *current_structure;

  if (ligand) {
    // look for second mfe root node! one of two mfe has a star (ligand bound structure).
    for (ii = 1; ii <= nlmin; ii++) {
      if (truemin[ii] == 0)
        continue;

      unsigned long root_gradmin = ii;
      while (lmin[root_gradmin].father != 0)
        root_gradmin = lmin[root_gradmin].father;

      current_lm        = &lmin[root_gradmin];
      current_structure = unpack_my_structure(current_lm->structure);
      if (current_structure[strlen(current_structure) - 1] == '*') {
        if (current_lm->energy < star_mfe) {
          star_mfe        = current_lm->energy;
          star_mfe_index  = root_gradmin;
        }
      } else {
        if (current_lm->energy < normal_mfe) {
          normal_mfe        = current_lm->energy;
          normal_mfe_index  = root_gradmin;
        }
      }

      free(current_structure);
    }
  }

  unsigned long mfe_comp_index = 0;

  for (ii = 1; ii <= nlmin; ii++) {
    if (truemin[ii] == 0)
      continue;

    unsigned long root_gradmin = ii;
    while (lmin[root_gradmin].father != 0)
      root_gradmin = lmin[root_gradmin].father;

    if ((!ligand && root_gradmin == 1) ||
        (ligand && ((root_gradmin == star_mfe_index) || (root_gradmin == normal_mfe_index))))
      //ii is connected to the mfe tree
      mfe_component_minima[mfe_comp_index++] = truemin[ii];
    else
      continue;
  }
  mfe_component_minima[mfe_comp_index] = 0; // zero terminated
  return mfe_component_minima;
}


void
ps_tree(loc_min       *Lmin,
        unsigned long *truemin)
{
  nodeT         *nodes;
  unsigned long i, ii;
  unsigned long nlmin;

  nlmin = Lmin[0].fathers_pool;

  if (max_print > truemin[0])
    max_print = truemin[0];

  nodes = (nodeT *)space(sizeof(nodeT) * (max_print + 1));
  for (i = 0, ii = 1; i < max_print && ii <= nlmin; ii++) {
    unsigned long s1, f;
    double        E_saddle;
    if ((s1 = truemin[ii]) == 0)
      continue;

    if (s1 > max_print)
      nrerror("inconsistency in ps_tree, aborting");

    E_saddle  = Lmin[ii].E_saddle;
    f         = Lmin[ii].father;
    if (f == 0)
      E_saddle = Lmin[0].E_saddle;         /* maximum energy */

    nodes[s1 - 1].father = (f == 0) ? ULONG_MAX : truemin[f] -1;

    nodes[s1 - 1].height        = Lmin[ii].energy;
    nodes[s1 - 1].saddle_height = E_saddle;


    if (print_labels) {
      char  *L;
      char  *s;
      s = unpack_my_structure(Lmin[ii].structure);
      if ((POV_size) && (Lmin[ii].global)) {
        L = (char *)space(sizeof(char) * (3 + strlen(s)));
        strcat(L, s);
        strcat(L, " *");
        nodes[s1 - 1].label = L;
      } else {
        nodes[s1 - 1].label = strdup(s);
      }

      free(s);
    } else {
      char *L = NULL, *s = NULL;
      L = (char *)space(sizeof(char) * 10);
      s = unpack_my_structure(Lmin[ii].structure);
      if (s[strlen(s) - 1] == '*') {
        (void)sprintf(L, "%ld", s1);
        nodes[s1 - 1].label = L;
      } else {
        free(L);
      }

      if ((POV_size) && (Lmin[ii].global)) {
        (void)sprintf(L, "%ld *", s1);
        nodes[s1 - 1].label = L;
      }

      free(s);
    }

    i++;
  }

  PS_tree_plot(nodes, max_print, "tree.ps");

  for (i = 0; i < (max_print); i++)
    if (nodes[i].label != NULL)
      free(nodes[i].label);

  free(nodes);
}


/*=========================*/
char *
get_taxon_label(int i)
{
  char *label;

  if (i == -1)
    return NULL;

  label = (char *)space(20);
  sprintf(label, "%d", i);

  return label;
}


/*===============================*/


static path_entry     *path;
static unsigned long  np, max_path;

static int
path_cmp(const void *a,
         const void *b)
{
  path_entry  *A, *B;
  int         d;

  A = (path_entry *)a;
  B = (path_entry *)b;
  if ((d = strcmp(A->key, B->key)) == 0)
    return A->num - B->num;
  else
    return d;
}


/*=======*/
path_entry *
backtrack_path(unsigned long  l1,
               unsigned long  l2,
               loc_min        *LM,
               unsigned long  *truemin,
               hash_table_t *hash_table)
{
  unsigned long n_lmin, i, ll1 = 0, ll2 = 0;
  char          *tag;

  lmin    = LM;
  n_lmin  = lmin[0].fathers_pool;
  for (i = 1; i <= n_lmin; i++) {
    if (truemin[i] == l1)
      ll1 = i;

    if (truemin[i] == l2)
      ll2 = i;
  }
  np        = 0;
  max_path  = 128;
  path      = (path_entry *)space(max_path * sizeof(path_entry));
  tag       = (char *)space(16);
  backtrack_path_rec(ll1, ll2, tag, hash_table);
  path[np].hp = NULL;
  qsort(path, np, sizeof(path_entry), path_cmp);
  free(tag);
  return path;
}


static void
backtrack_path_rec(unsigned long  l1,
                   unsigned long  l2,
                   const char     *tag,
                   hash_table_t* hash_table)
{
  hash_entry    h, *l1dir = NULL, *l2dir = NULL;
  int           dir = 1;
  unsigned long child, father, maxsaddle;

  /* if left==1 left points toward l2 else toward l1 */
  if (l1 > l2) {
    dir = -1;
    {
      unsigned long t;
      t   = l1;
      l1  = l2;
      l2  = t;
    }
  }

  child   = l2;
  father  = l1;

  maxsaddle = child;
  /* find saddle connecting l1 and l2 */
  while (lmin[child].father != father) {
    if (lmin[child].father == 0) {
      fprintf(stderr, "ERROR in backtrack_path(): ");
      fprintf(stderr, "No saddle between lmin %ld and lmin %ld\n", l2, l1);
      exit(1);
    }

    child = lmin[child].father;
    if (child < father) {
      unsigned long t;
      t       = child;
      child   = father;
      father  = t;
    }

    if (lmin[child].E_saddle > lmin[maxsaddle].E_saddle)
      maxsaddle = child;

    /* fprintf(stderr,"f:>%d< c:>%d< %d\n", father, child, maxsaddle); */
  }
  /* found the saddle point, maxsaddle, connecting l1 and l2 */
  h.structure = lmin[maxsaddle].saddle;
  //path[np].hp = lookup_hash(&h);
  path[np].hp = (hash_entry *)ht_get(*hash_table, &h);
  strcpy(path[np].key, tag);
  strcat(path[np].key, "M");
  np++;
  if (np + 2 >= max_path) {
    max_path  *= 2;
    path      = (path_entry *)xrealloc(path, max_path * sizeof(path_entry));
  }

  /* which direction from saddle to l2, l1 ? */
  for (child = l2; child > 0; child = lmin[child].father) {
    if (child == maxsaddle) {
      l2dir = lmin[maxsaddle].left;
      l1dir = lmin[maxsaddle].right;
      break;
    }

    if (child == lmin[maxsaddle].father) {
      l2dir = lmin[maxsaddle].right;
      l1dir = lmin[maxsaddle].left;
      break;
    }
  }

  /* branch to l2,  else saddle==l2 and we're done */
  if (l2dir)
    walk_limb(l2dir, l2, -dir, tag, hash_table);

  /* branch to l1 (to father) */
  if (l1dir) /* else saddle==l1 and we're done */
    walk_limb(l1dir, l1, dir, tag, hash_table);
}


/*=======================================================================*/
static void
walk_limb(hash_entry    *hp,
          unsigned long LM,
          int           inc,
          const char    *tag,
          hash_table_t *hash_table)
{
  char          *tmp;
  int num = 0;
  hash_entry    *htmp;

  tmp = (char *)space(strlen(tag) + 4);

  strcpy(tmp, tag);
  strcat(tmp, (inc > 0) ? "R" : "LZ");
  /* walk down until u hit a local minimum */
  for (htmp = hp; htmp->down != NULL; htmp = htmp->down, num += inc, np++) {
    if (np + 2 >= max_path) {
      max_path  *= 2;
      path      = (path_entry *)xrealloc(path, max_path * sizeof(path_entry));
    }

    path[np].hp = htmp;
    strcpy(path[np].key, tmp);
    path[np].num = num;
  }
  /* store local minimum (but only once) */
  if (htmp->basin == LM) {
    path[np].hp = htmp;
    strcpy(path[np].key, tmp);
    path[np++].num = num;
    if (np + 2 >= max_path) {
      max_path  *= 2;
      path      = (path_entry *)xrealloc(path, max_path * sizeof(path_entry));
    }
  }

  if (inc < 0)
    tmp[strlen(tmp) - 1] = '\0';

  /* fprintf(stderr, "walk towards %d lands in %d\n", LM, htmp->basin); */

  /* wrong local minimum start cruising again */
  if (htmp->basin != LM) {
    if (inc == -1)
      backtrack_path_rec(htmp->basin, LM, tmp, hash_table);
    else
      backtrack_path_rec(LM, htmp->basin, tmp, hash_table);
  }

  free(tmp);
}


void
print_path(FILE           *PATH,
           path_entry     *path,
           unsigned long  *tm,
           unsigned long  *mfe_component_true_min_indices)
{
  unsigned long i, n, grad_min_index, mfe_component_index;
  int           found_gradmin;

  for (i = 0; path[i].hp; i++) {
    char c[32] = {
      0
    }, *struc;
    if (path[i].hp->down == NULL && tm[path[i].hp->basin] != 0) {
      grad_min_index = tm[path[i].hp->basin];

      if (mfe_component_true_min_indices == NULL) {
        sprintf(c, "L%04ld", grad_min_index);
      } else {
        /* find gradient minimum in mfe component index list*/
        n             = 0;
        found_gradmin = 0;
        while (mfe_component_true_min_indices[n] != 0) {
          if (grad_min_index == mfe_component_true_min_indices[n]) {
            mfe_component_index = n + 1;
            found_gradmin       = 1;
            sprintf(c, "L%04ld", mfe_component_index);
            break;
          }

          n++;
        }
        if (found_gradmin == 0) {
          fprintf(stderr,
                  "Error: the minimum on the path is not in the connected component of the mfe structure! The path is not complete!\n");
          break;
        }
      }
    } else
    if (path[i].key[strlen(path[i].key) - 1] == 'M') {
      c[0] = 'S';
    } else {
      c[0] = 'I';
    }

    struc = unpack_my_structure(path[i].hp->structure);
    fprintf(PATH, "%s (%6.2f) %-5s\n", struc, path[i].hp->energy, c);
    free(struc);
  }
}


static void
merge_components(unsigned long  c1,
                 unsigned long  c2)
{
  if (comp[c1].size < comp[c2].size) {
    unsigned long cc;
    cc  = c1;
    c1  = c2;
    c2  = cc;
  }

  comp[c1].size += comp[c2].size;
  truecomp[c2]  = c1;
  set_merge(comp[c1].basins, comp[c2].basins);
}


static int
comp_comps(const void *A,
           const void *B)
{
  struct comp   *a, *b;
  unsigned long i = 0, ba, bb;

  a = (struct comp *)A;
  b = (struct comp *)B;
  for (i = 0; i < a->basins->num_elem && i < b->basins->num_elem; i++) {
    ba = a->basins->data[i].basin;
    bb = b->basins->data[i].basin;
    if (ba != bb){
      if (ba > bb)
        return 1;
      else
        return -1;
    }
  }
  return (i == a->basins->num_elem) ? -1 : 1;
}


map_struc
get_mapstruc(char           *p,
             loc_min        *LM,
             unsigned long  *tm,
             hash_table_t* hash_table)
{
  hash_entry    *hp, h;
  char          *pp;
  unsigned long min, gradmin;
  map_struc     ms;

  ms.structure  = NULL;
  pp            = pack_my_structure(p);
  h.structure   = pp;
  //hp            = lookup_hash(&h);
  hp= (hash_entry *)ht_get(*hash_table, &h);

  if (hp == NULL) {
    fprintf(stderr, "get_mapstruc: structure not in hash\n");
    free(pp);
    return ms;
  }

  min = hp->basin;
  while (tm[min] == 0)
    min = LM[min].father;
  gradmin = hp->GradientBasin;
  while (tm[gradmin] == 0)
    gradmin = LM[gradmin].father;

  if (gradmin == 0) {
    fprintf(stderr, "get_mapstruc: gradient minimum not yet assigned\n");
    free(pp);
    return ms;
  }

  ms.structure          = unpack_my_structure(LM[gradmin].structure);
  ms.n                  = hp->n;
  ms.min                = min;
  ms.truemin            = tm[min];
  ms.gradmin            = gradmin;
  ms.truegradmin        = tm[gradmin];
  ms.energy             = hp->energy;
  ms.min_energy         = LM[min].energy;
  ms.truemin_energy     = LM[tm[min]].energy;
  ms.gradmin_energy     = LM[gradmin].energy;
  ms.truegradmin_energy = LM[tm[gradmin]].energy;
  free(pp);
  return ms;
}


void
print_rates(unsigned long       n,
            barrier_options     *opt,
            barriers_rates_type rate_files)
{
  unsigned long i, j;
  FILE          *OUT;

  if (rate_files & Barriers_binary_rates) {
    FILE    *BINOUT;
    char    *binfile = opt->binary_rates_file;
    double  tmprate;
    BINOUT = fopen(binfile, "w");
    if (!BINOUT) {
      fprintf(stderr, "could not open file pointer 4 binary outfile\n");
      exit(101);
    }

    if (n > INT_MAX) {
      fprintf(stderr, "Error: more local minima than 32bit integer values!");
      // we need a new binary format in order to support this.
      exit(EXIT_FAILURE);
    }

    /* first write dim to file */
    fwrite(&n, sizeof(int), 1, BINOUT);
    for (i = 1; i <= n; i++)
      for (j = 1; j <= n; j++) {
        tmprate = rate[j][i];
        fwrite(&tmprate, sizeof(double), 1, BINOUT);
      }
    fprintf(stderr, "rate matrix written to binfile\n");
    fclose(BINOUT);
  }

  OUT = fopen(opt->text_rates_file, "w");
  if (!OUT) {
    fprintf(stderr, "could not open rates file %s for output\n", opt->text_rates_file);
    return;
  }

  for (i = 1; i <= n; i++) {
    for (j = 1; j <= n; j++)
      fprintf(OUT, "%10.4g ", rate[i][j]);
    fprintf(OUT, "\n");
    free(rate[i]);
  }
  free(rate);
  fclose(OUT);
}


void
compute_rates(unsigned long *truemin,
              char          *sequence,
              hash_table_t* hash_table)
{
  unsigned long i, j, ii, r, gb, gradmin, n, rc, *realnr = NULL;
  char          *p, *pp, *structure, newsub[10] = "new.sub", mr[15] = "microrates.out";
  hash_entry    *hpr, h, *hp;
  double        Zi;
  FILE          *NEWSUB = NULL, *MR = NULL;

  if (ligand == 1)
    move_it = RNA_move_it_rates;

  n     = truemin[0];
  rate  = (double **)space((n + 1) * sizeof(double *));
  dr    = (double *)space((n + 1) * sizeof(double));
  for (i = 1; i <= n; i++)
    rate[i] = (double *)space((n + 1) * sizeof(double));
  if (do_microrates) {
    realnr  = (unsigned long *)space((Read_lines + 1) * sizeof(unsigned long));
    MR      = fopen(mr, "w");
    NEWSUB  = fopen(newsub, "w");
    fprintf(NEWSUB, "%s %6.2f\n", sequence, 100 * mfe);
    fflush(NEWSUB);
    fprintf(MR, ">%ld states\n", Read_lines);
  }

  for (rc = 1, r = 0; r < Read_lines; r++) {
    unsigned long b;
    hpr     = &hpool[r];
    Zi      = exp((mfe - hpr->energy) / kT);
    gradmin = hpr->GradientBasin;
    while (truemin[gradmin] == 0) {
      gradmin = lmin[gradmin].father;
      if (gradmin == 0)
        break; // break if we are at the root!
    }
    gradmin = truemin[gradmin];
    if (gradmin > n)
      continue;

    for (b = hpr->basin; b > 1; b = lmin[b].father);
    structure = unpack_my_structure(hpr->structure);
    move_it(structure);       /* generate all neighbors of configuration */

    for (i = 0; i <= n; i++)
      dr[i] = 0;
    while ((p = pop())) {
      int double_move = 0;
      if (IS_RNA) {
        if (p[strlen(p) - 1] == 'D') {
          double_move       = 1;
          p[strlen(p) - 1]  = '\0';
        } else {
          double_move = 0;
        }
      }

      pp          = pack_my_structure(p);
      h.structure = pp;
      /* check whether we've seen the structure before */
      //if ((hp = lookup_hash(&h))) {
      if ((hp= (hash_entry *)ht_get(*hash_table, &h))) {
        if (hp->n <= r) {
          gb = hp->GradientBasin;
          while (truemin[gb] == 0)
            gb = lmin[gb].father;
          gb = truemin[gb];
          if (gb <= n)
            dr[gb] += (double_move) ? (noLP_rate * Zi) : Zi;

          if (do_microrates && b && (realnr)) {
            double rate, dg;
            dg    = hpr->energy - hp->energy;
            rate  = exp(-dg / kT);
            fprintf(MR, "%10ld %8ld %15.12f 1\n", rc, realnr[hp->n], rate);
          }
        }
      }

      free(pp);
    }
    if (do_microrates && b) {
      fprintf(NEWSUB, "%s %6.2f %li %li\n", structure, hpr->energy, gradmin, hpr->basin);
      fflush(NEWSUB);
      realnr[hpr->n] = rc++;
    }

    for (i = 1; i <= n; i++) {
      rate[i][gradmin]  += dr[i];
      rate[gradmin][i]  += dr[i];
    }
    free(structure);
    reset_stapel();
  }

  fprintf(stderr, "done with 2nd pass\n");
  free(dr);

  for (i = ii = 1; i <= n; i++, ii++) {
    while (truemin[ii] != i)
      ii++;
    for (j = 1; j <= n; j++)
      rate[i][j] /= lmin[ii].Zg;
  }
  if (do_microrates) {
    free(realnr);
    fclose(NEWSUB);
    fclose(MR);
  }

  if (free_move_it)
    free_move_it();

  free_stapel();
}


void
free_rates(unsigned long length_rates)
{
  unsigned long i;

  for (i = 0; i <= length_rates; i++)
    free(rate[i]);
  free(rate);
}


void
print_rates_of_mfe_component(unsigned long        *mfe_component_true_min_indices,
                             barrier_options      *opt,
                             barriers_rates_type  rate_files)
{
  unsigned long n, i, j, ii, jj;
  FILE          *OUT;

  if (mfe_component_true_min_indices != NULL) {
    n = 0;
    while (mfe_component_true_min_indices[n] != 0)
      n++;
  } else {
    fprintf(stderr, "Error: could not print rates because mfe component is NULL!");
    return;
  }

  if (rate_files & Barriers_binary_rates) {
    FILE    *BINOUT;
    char    *binfile = opt->binary_rates_file;
    double  tmprate;
    BINOUT = fopen(binfile, "w");
    if (!BINOUT) {
      fprintf(stderr, "could not open file pointer 4 binary outfile\n");
      exit(101);
    }

    if (n > INT_MAX) {
      fprintf(stderr, "Error: more local minima than 32bit integer values!");
      // we need a new binary format in order to support this.
      exit(EXIT_FAILURE);
    }

    /* first write dim to file */
    fwrite(&n, sizeof(int), 1, BINOUT);
    for (i = 0; i < n; i++) {
      ii = mfe_component_true_min_indices[i];
      for (j = 0; j < n; j++) {
        jj      = mfe_component_true_min_indices[j];
        tmprate = rate[jj][ii];
        fwrite(&tmprate, sizeof(double), 1, BINOUT);
      }
    }
    fprintf(stderr, "rate matrix written to binfile\n");
    fclose(BINOUT);
  }

  OUT = fopen(opt->text_rates_file, "w");
  if (!OUT) {
    fprintf(stderr, "could not open rates file %s for output\n", opt->text_rates_file);
    return;
  }

  for (i = 0; i < n; i++) {
    ii = mfe_component_true_min_indices[i];
    for (j = 0; j < n; j++) {
      jj = mfe_component_true_min_indices[j];
      fprintf(OUT, "%10.4g ", rate[ii][jj]);
    }
    fprintf(OUT, "\n");
  }
  fclose(OUT);
}


void
ps_tree_mfe_component(loc_min       *Lmin,
                      unsigned long *truemin,
                      unsigned long *mfe_component_true_min_indices)
{
  nodeT         *nodes;
  unsigned long i, ii;
  unsigned long nlmin;

  nlmin = Lmin[0].fathers_pool;

  //if (max_print > truemin[0])
  //   max_print = truemin[0];
  unsigned long mfe_comp_max;
  for (mfe_comp_max = 0; mfe_component_true_min_indices[mfe_comp_max] != 0; mfe_comp_max++);
  unsigned long max_print = mfe_comp_max;

  nodes = (nodeT *)space(sizeof(nodeT) * (max_print + 1));
  for (i = 0, ii = 1; i < max_print && ii <= nlmin; ii++) {
    unsigned long s1, f;
    double        E_saddle;
    if ((s1 = truemin[ii]) == 0)
      continue;

    unsigned long mfe_comp_index = 0;
    for (mfe_comp_index = 0; mfe_component_true_min_indices[mfe_comp_index] != 0; mfe_comp_index++)
      if (s1 == mfe_component_true_min_indices[mfe_comp_index])
        break;

    if (mfe_comp_index == mfe_comp_max)
      continue;

    if (i > max_print)
      nrerror("inconsistency in ps_tree, aborting");

    E_saddle  = Lmin[ii].E_saddle;
    f         = Lmin[ii].father;
    if (f == 0)
      E_saddle = Lmin[0].E_saddle;         /* maximum energy */

    unsigned long mfe_comp_father_index = 0;
    for (mfe_comp_father_index = 0; mfe_component_true_min_indices[mfe_comp_father_index] != 0;
         mfe_comp_father_index++)
      if (truemin[f] == mfe_component_true_min_indices[mfe_comp_father_index])
        break;

    if (f != 0 && mfe_comp_father_index >= mfe_comp_max)
      continue;                                                           // do not include it if father is not the root and not in the mfe component.

    nodes[mfe_comp_index].father = (f == 0) ? ULONG_MAX : mfe_comp_father_index; //(f == 0) ? -1 : truemin[f] - 1;

    nodes[mfe_comp_index].height        = Lmin[ii].energy;
    nodes[mfe_comp_index].saddle_height = E_saddle;

    if (print_labels) {
      char  *L;
      char  *s;
      s = unpack_my_structure(Lmin[ii].structure);
      if ((POV_size) && (Lmin[ii].global)) {
        L = (char *)space(sizeof(char) * (3 + strlen(s)));
        strcat(L, s);
        strcat(L, " *");
        nodes[mfe_comp_index].label = L;
      } else {
        nodes[mfe_comp_index].label = strdup(s);
      }

      free(s);
    } else {
      char *L = NULL, *s = NULL;
      L = (char *)space(sizeof(char) * 10);
      s = unpack_my_structure(Lmin[ii].structure);
      if (s[strlen(s) - 1] == '*') {
        (void)sprintf(L, "%ld", s1);
        nodes[mfe_comp_index].label = L;
      } else {
        free(L);
      }

      if ((POV_size) && (Lmin[ii].global)) {
        (void)sprintf(L, "%ld *", s1);
        nodes[mfe_comp_index].label = L;
      }

      free(s);
    }

    i++;
  }

  PS_tree_plot(nodes, max_print, "tree.ps");

  for (i = 0; i < (max_print); i++)
    if (nodes[i].label != NULL)
      free(nodes[i].label);

  free(nodes);
}
