/* Last changed Time-stamp: <2001-04-08 00:52:46 ivo> */
/* barriers.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <math.h>
#include <limits.h>
#include "ringlist.h"
#include "stapel.h"
#include "utils.h"
#include "hash_util.h"
#include "barrier_types.h"
#include "treeplot.h"

/* Tons of static arrays in this one! */
static char UNUSED rcsid[] = "$Id: barriers.c,v 1.2 2001/04/09 08:02:00 ivo Exp $";
#ifdef __GNUC__BLAH
#define XLL unsigned long long
#define MAXIMUM 18446744073709551615ULL
#else
#define XLL unsigned int
#define MAXIMUM ULONG_MAX
#endif

static char *form;         /* array for configuration */ 
static loc_min *lmin;      /* array for local minima */

/* "global" variables */

static int n_lmin;
static unsigned int max_lmin;
static int max_print;
static int n_saddle;
static double minh=0.0000001;
static double energy;
    /* energy of last read structure (for check_neighbors) */
static double mfe;           /* used for scaling Z */

static void (*move_it)(char *);
static void (*free_move_it)(void) = NULL;
static char *(*pack_my_structure)(const char *) ;
static char *(*unpack_my_structure)(const char *) ;

static double kT= -1;
extern unsigned long collisions;

/* global switches */  /* defaults changed */
static int print_saddles = 1;
static int bsize = 1;
static int shut_up = 0;
static int verbose = 0;
static int max_print = 0;
static int IS_RNA = 0;

/* private functions */
static void walk_limb (hash_entry *hp, int LM, int inc, const char *tag);
static void backtrack_path_rec (int l1, int l2, const char *tag);
static int *make_sorted_index(int *truemin);
static void Sorry(char *GRAPH);
static int  read_data(FILE *FP, double *energy,char *strucb, int len);

/* public functiones */
int      *make_truemin(loc_min *Lmin);
loc_min  *barriers(barrier_options opt);

static int  compare(const void *a, const void *b);
void check_neighbors(void);
void merge_basins(void);
void print_results(loc_min *L, int *tm);
void ps_tree(loc_min *Lmin, int *truemin);
     
/* ----------------------------------------------------------- */

void set_barrier_options(barrier_options opt) {
  print_saddles = opt.print_saddles;
  bsize = opt.bsize;
  shut_up = opt.want_quiet;
  max_print = opt.max_print;
  minh = opt.minh;
  verbose = opt.want_verbose;
  switch(opt.GRAPH[0]) {
  case 'R' :    /* RNA secondary Structures */
    if (strncmp(opt.GRAPH, "RNA", 3)==0) {
      int nolp=0, shift=1;
      IS_RNA=1;
      if (opt.kT<=-300) opt.kT=37;
      kT = 0.00198717*(273.15+opt.kT);   /* kT at 37C in kcal/mol */
      move_it = RNA_move_it;
      free_move_it = RNA_free_rl;
      pack_my_structure = pack_structure;
      unpack_my_structure = unpack_structure;
      if (strstr(opt.GRAPH,   "noLP")) nolp=1;
      if (strstr(opt.MOVESET, "noShift")) shift=0;
      RNA_init(opt.seq, shift, nolp);
      if (verbose) 
	fprintf(stderr, "Graph is RNA with noLP=%d, Shift=%d\n", nolp, shift);
    } else Sorry(opt.GRAPH);
    break;
  case 'Q' :    /* Haming graphs */
    if (opt.GRAPH[1]=='2') {   /* binary +- alphabet */
      move_it = SPIN_move_it;
      pack_my_structure = strdup;
      unpack_my_structure = strdup;
      if (verbose) 
	fprintf(stderr, "Graph is Q2\n");
    }
    else Sorry(opt.GRAPH);
    break;
  case 'P' :    /* Permutations */
    switch(*opt.MOVESET) {
    case 'R' :
      move_it = Reversal_move_it; break;
    case 'C' :
      move_it = CTranspos_move_it; break;
    case 'T':
    default:
      
      move_it = Transpos_move_it;
    }
    pack_my_structure = strdup;
    unpack_my_structure = strdup;
    if (verbose) 
      fprintf(stderr, "Graph is Permutations with moveset %c\n",
	      *opt.MOVESET ? *opt.MOVESET : 'T');
    break;
  case 'T' :    /* Phylogenetic Trees */
    move_it = NNI_move_it;
    pack_my_structure = strdup;
    unpack_my_structure = strdup;
    if (verbose) 
      fprintf(stderr, "Graph is Trees with NNI moves\n");
    break;   
  default :
    Sorry(opt.GRAPH);
  }
  if (kT<0) {
    if (opt.kT<=-300) kT=1;
    else kT=opt.kT;
  }
}

static void Sorry(char *GRAPH) {
  fprintf(stderr,"Graph \"%s\" is not implemented\n",GRAPH);
  exit(-2);
}

loc_min *barriers(barrier_options opt) {
  int length, read =0;

  set_barrier_options(opt);
  length = strlen(opt.seq);
  max_lmin = 1023;
  lmin = (loc_min *) space((max_lmin + 1) * sizeof(loc_min));
  n_lmin = 0;

  form = (char *) space((length+1)*sizeof(char));
  
  ini_stapel(length);
  
  while (read_data(opt.INFILE, &energy,form,length)) {
    if (read==0) mfe=energy;
    read++;   
    move_it(form);       /* generate all neighbor of configuration */
    check_neighbors();   /* flood the energy landscape */
    reset_stapel();
    if (n_saddle+1 == max_print)
      break;  /* we've found all we want to know */
  }
  
  if(!shut_up) fprintf(stderr,
		       "read %d structures, to find %d saddles\n",
		       read, n_saddle);
  
  if (max_print == 0 || max_print > n_lmin)
    max_print = n_lmin;
  
  lmin[0].fathers_pool = n_lmin;   /* store size here; pfs 03 2001 */
  lmin[0].E_saddle = energy + 0.001;
  lmin[0].energy = lmin[1].energy;

  if (free_move_it) 
    free_move_it();
  free_stapel();
  free(form);
  fflush(stdout);
  if(!shut_up) fprintf(stderr, "%ld hash table collisions\n", collisions);

  return lmin;
}

int *make_truemin(loc_min *Lmin) {
  int *truemin;
  int nlmin;
  nlmin = Lmin[0].fathers_pool;
  truemin = (int *) space((nlmin+1)*sizeof(int));
  /* truemin[0] = nlmin; */
  {
    int i,ii;
    for (ii=i=1; (i<=max_print)&&(ii<=n_lmin); ii++) {
      if (!lmin[ii].father) lmin[ii].E_saddle = energy + 0.000001;
      if (lmin[ii].E_saddle - lmin[ii].energy > minh) {
	truemin[ii]=i++;
      }
    }
  }
  return truemin;
}



/*=============================================================*/
static int read_data(FILE *FP, double *energy,char *strucb, int len) {

  static char struc[501];
  int r;

  r = fscanf(FP, "%500s %lf", struc, energy);
  if (r==EOF) return 0;
  if (r!=2) {
    fprintf(stderr, "Error in input file\n");
    exit(123);
  }

  if(strlen(struc) != len) {
    fprintf(stderr,"read_data():\n%s\n unequal length !!\n", struc);
    exit (1);
  }
  strcpy(strucb, struc);

  return (1);
}

typedef struct {
  int basin;
  hash_entry *hp;
} basinT;

/*=====================================*/
static int compare(const void *a, const void *b) {
  int A, B;
  A = ((basinT *)a)->basin; B = ((basinT *)b)->basin;
  return (A - B);
}


/*======================*/
void check_neighbors(void)
{
  char *p, *pp, *pform;
  int nb, i_lmin, basin, obasin;
  static int false_lmin=0;
  hash_entry *hp, h, *down=NULL;
  basinT basins[1000];
    
  float minenergia;         /* for Gradient Basins */
  int   gradmin=0;            /* for Gradient Basins */
  
  nb = 0; obasin = -1;
  minenergia = 100000000.0; /* for Gradient Basins */
  
  /* foreach neighbor structure of configuration "Structure" */
  while ((p = pop()))
    {
      pp = pack_my_structure(p); 
      h.structure = pp;

      /* check whether we've seen the structure before */
      if ((hp = lookup_hash(&h)))
	{
	  /* because we've seen this structure before, it */
	  /* already belongs to the basin of attraction */
	  /* of a local minimum */
	  basin = hp->basin;

	  if ( hp->energy < minenergia ) { /* for Gradient Basins */
	    minenergia = hp->energy;       /* for Gradient Basins */
	    gradmin = hp->GradientBasin;   /* for Gradient Basins */
	    down = hp;
	  }                                /* for Gradient Basins */
 
	  /* the basin of attraction of this local minimum may have been */
	  /* merged with the basin of attraction of an energetically */
	  /* "deeper" local minimum in a previous step */
	  /* go and find this "deeper" local minimum! */
	  while (lmin[basin].father) basin=lmin[basin].father;
	  
	  /* put the "deepest" local minimum into the basins-list */
	  if (basin != obasin) {
	    basins[nb].hp = hp;
	    basins[nb++].basin = basin;
	  }
	  obasin = basin;
	}
      free(pp);
    }

  pform = pack_my_structure(form);

  if (nb > 1) {
    /* Structure belongs to at least 2 basins of attraction */
    /* i.e. it's a saddle */
    int i;
    int pool = 0;
    double Z=0;

    /* sort the basins-list in ascending order */
    qsort(basins, nb, sizeof(basinT), compare);
    
    /* merge all basins of attraction Structure belongs to */
    for (i = 1; i < nb; i++)
      {
	int ii;
	ii = basins[i].basin;
	if (ii == basins[i-1].basin) /* been there */
	  continue;
	/* going to merge basin[i] with basin[0] */
	if ((!max_print) || (ii<=max_print+false_lmin)) { 
	  /* found the saddle for a basin we're gonna print */
	  if (energy-lmin[ii].energy>minh) n_saddle++;
	  else false_lmin++;
	}
	
	if (lmin[ii].father != 0) fprintf(stderr, "This shouldn't happen\n");
	lmin[ii].father = basins[0].basin;
	lmin[ii].saddle = pform;
	lmin[ii].E_saddle = energy;
	lmin[ii].left = basins[i].hp;
        lmin[ii].right = basins[0].hp;
	if (bsize)
	  {
	    lmin[ii].fathers_pool = lmin[basins[0].basin].my_pool;
	    pool += lmin[ii].my_pool;
	    Z += lmin[ii].Z; 
	  }
      }
    if (bsize) {
      lmin[basins[0].basin].my_pool += pool;
      lmin[basins[0].basin].Z += Z;
    }
  }
  
  if (nb>0)
    i_lmin = basins[0].basin;

  else {
    /* Structure is a "new" local minimum */
    
    i_lmin = n_lmin+1;
    n_lmin = i_lmin;

    gradmin = n_lmin;        /* for Gradient Basins */
    down = NULL;
    /* need to allocate more space for the lmin-list */
    if (n_lmin > max_lmin) {
      fprintf(stderr, "increasing lmin array\n");
      lmin = (loc_min *) xrealloc(lmin, (max_lmin*2+1)*sizeof(loc_min));
      memset(lmin + max_lmin +1, 0, max_lmin);
      max_lmin *= 2;
    }
    
    /* store configuration "Structure" in lmin-list */
    lmin[i_lmin].father = 0;
    lmin[i_lmin].structure = pform;
    lmin[i_lmin].energy = energy;
    lmin[i_lmin].my_GradPool = 0;
    lmin[i_lmin].Z = 0;
    lmin[i_lmin].Zg = 0;
  }

  /* store configuration "Structure" in hash table */
  hp = (hash_entry *) space(sizeof(hash_entry));
  hp->structure = pform;
  hp->energy = energy;
  hp->basin = i_lmin;
  hp->GradientBasin = gradmin;    /* for Gradient Basins */
  hp->down = down;
  lmin[i_lmin].my_pool++;
  lmin[i_lmin].Z += exp(-(energy-mfe)/kT);
  lmin[gradmin].my_GradPool++;
  lmin[gradmin].Zg += exp(-(energy-mfe)/kT);
  if (energy<mfe) fprintf(stderr, "%f %f shouldn't happen!\n", energy, mfe);
  if (write_hash(hp))
    nrerror("duplicate structure");
}

/*====================*/
void print_results(loc_min *Lmin, int *truemin)
{
  int i,ii,j;
  char *struc;
  char *format;
  
  if (IS_RNA)
    format = "%4d %s %6.2f %4d %6.2f";
  else
    format = "%4d %s %13.5f %4d %13.5f";
      
  n_lmin = Lmin[0].fathers_pool;

  /* printf("     %s\n", farbe); */
  for (i = 1; i <= n_lmin; i++)
    {
      if ((ii = truemin[i])==0) continue;

      struc = unpack_my_structure(Lmin[i].structure);
      printf(format, ii, struc, Lmin[i].energy,
	     truemin[Lmin[i].father], Lmin[i].E_saddle - Lmin[i].energy);
      free(struc);
      
      if (print_saddles) { 
	if (Lmin[i].saddle)  {
	  struc = unpack_my_structure(Lmin[i].saddle);
	  printf(" %s", struc);
	  free(struc);
	}
        else {
	  printf(" ");
	  for (j=0;j<strlen(struc);j++) { printf("~"); }
	}
	 
      }
      if (bsize) 
	printf (" %12ld %8ld %7.3f %8ld %7.3f",
		Lmin[i].my_pool, Lmin[i].fathers_pool, mfe -kT*log(lmin[i].Z),
		Lmin[i].my_GradPool, mfe -kT*log(lmin[i].Zg));
      printf("\n");
    }
}

typedef struct {
  int   set1;
  int   set2;
  float distance;
  float distance2;
} Union;

void ps_tree(loc_min *Lmin, int *truemin)
{
  nodeT *nodes;
  int i,ii;
  int nlmin;
  
  nlmin = Lmin[0].fathers_pool;
  
  if (max_print>nlmin) max_print=nlmin;
  
  nodes = (nodeT *) space(sizeof(nodeT)*(max_print+1));
  for (i=0,ii=1; i<max_print && ii<=nlmin; ii++)
    {
      register int s1, f;
      double E_saddle;
      if ((s1=truemin[ii])==0) continue;
      if (s1>max_print) 
	nrerror("inconsistency in ps_tree, aborting");
      E_saddle = Lmin[ii].E_saddle;
      f = Lmin[ii].father; 
      if (f==0) {
	E_saddle = Lmin[0].E_saddle; /* maximum energy */
	f=1;                         /* join with mfe  */
      }
      nodes[s1-1].father = truemin[f]-1;
      nodes[s1-1].height = Lmin[ii].energy;
      nodes[s1-1].saddle_height = E_saddle;
      i++;
    }
  PS_tree_plot(nodes, max_print, "tree.ps");
  free(nodes);
}

/*=========================*/
char * get_taxon_label(int i)
{
  char *label;

  if (i == -1)
    return (NULL);
  label = (char *) space(20);
  sprintf(label, "%d", i);

  return (label);
}

/*=========================================*/
static int indx_comp(const void *i1, const void *i2)
{
  int i,j;

  i = *((int *)i1);
  j = *((int *)i2);
  if (lmin[i].E_saddle == lmin[j].E_saddle)
    return (i-j);
  return (lmin[i].E_saddle < lmin[j].E_saddle) ? (-1) : (1);
}

/*===============================*/
static int *make_sorted_index(int *truemin)
{
  int i, ii, *index;

  index = (int *) space((max_print+1)*sizeof(int));

  /* include only up to max_print local minima */
  for (i = 0, ii=1; i < max_print; ii++)
    if (truemin[ii])
      index[i++] = ii;
  
  /* sort local minima by saddle-point-heights */
  /* starting with the 1st excited state */
  qsort(index+1, max_print-1, sizeof(int), indx_comp);

  return (index);
}

static path_entry *path;
static int np, max_path=128;

static int path_cmp(const void *a, const void *b) {
  path_entry *A, *B; int d;
  A = (path_entry *) a;
  B = (path_entry *) b;
  if ((d=strcmp(A->key, B->key))==0) return (A->num-B->num);
  else return (d);
}

/*=======*/
path_entry *backtrack_path(int l1, int l2, loc_min *LM, int *truemin) {
  int n_lmin, i, ll1, ll2;
  char *tag;
  lmin = LM;
  n_lmin = lmin[0].fathers_pool;
  for (i=1; i<=n_lmin; i++) {
    if (truemin[i]==l1) ll1=i;
    if (truemin[i]==l2) ll2=i;
  }
  path = (path_entry *) space(max_path*sizeof(path_entry));
  tag = (char *) space(16);
  backtrack_path_rec(ll1, ll2, tag);  path[np].hp = NULL;
  qsort(path, np, sizeof(path_entry), path_cmp);

  return(path);
}

static void backtrack_path_rec (int l1, int l2, const char *tag)
{
  hash_entry h;
  int dir=1, left=1, swap=0, child, father, maxsaddle;
  /* if left==1 left points toward l2 else toward l1 */
  if (l1>l2) {
    dir = -1;
    {int t; t=l1; l1=l2; l2=t;}
  }
  child  = l2; father = l1;
  maxsaddle = child;
  /* find saddle connecting l1 and l2 */
  while (lmin[child].father != father) { 
    int tmp;
    if (lmin[child].father == 0){
      fprintf(stderr, "ERROR in backtrack_path(): ");
      fprintf(stderr,"No saddle between lmin %d and lmin %d\n", l2, l1);
      exit (1);
    }
    child = lmin[child].father;
    if (child<father) {tmp = child; child = father; father = tmp; swap= !swap;}
    if (lmin[child].E_saddle > lmin[maxsaddle].E_saddle) {
      maxsaddle = child;
      if (swap) left= -left;
    }
  }
  h.structure = pack_my_structure(lmin[maxsaddle].saddle);
  path[np].hp = lookup_hash(&h); free(h.structure);
  strcpy(path[np].key,tag); strcat(path[np].key, "M");
  np++;

  if (left>0) {
    /* branch to l2 */
    walk_limb (lmin[maxsaddle].left, l2, -dir, tag);
    /* reconstruct right branch (to father) */
    walk_limb (lmin[maxsaddle].right, l1, dir, tag);
  } else {
    /* branch to l2 */
    walk_limb (lmin[maxsaddle].right, l2, -dir, tag);
    /* reconstruct right branch (to father) */
    walk_limb (lmin[maxsaddle].left, l1, dir, tag);
  }    
}

/*=======================================================================*/
static void walk_limb (hash_entry *hp, int LM, int inc, const char *tag)
{
  char *tmp; int num=0;
  hash_entry *htmp;

  tmp = (char *) space(strlen(tag)+4);;
  strcpy(tmp, tag);
  strcat(tmp, (inc>0) ? "R" : "LZ");
  /* walk down until u hit a local minimum */
  for (htmp = hp; htmp->down != NULL; htmp = htmp->down, num += inc, np++) {
    if (np+2>=max_path) {
      max_path *= 2;
      path = (path_entry *) xrealloc(path, max_path*sizeof(path_entry));
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
  }

  if (inc<0) tmp[strlen(tmp)-1] = '\0';
  /* wrong local minimum start cruising again */
  if (htmp->basin != LM) {
    if (inc == -1)
      backtrack_path_rec (htmp->basin, LM, tmp);
    else
      backtrack_path_rec (LM, htmp->basin, tmp);
  }
}
