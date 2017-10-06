/* Last changed Time-stamp: <2017-10-06 16:27:13 mtw> */
/* barriers.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <stdbool.h>
#include "ringlist.h"
#include "stapel.h"
#include "utils.h"
#include "hash_util.h"
#include "barrier_types.h"
#include "compress.h"
#include "treeplot.h"
#include "simple_set.h"
#if HAVE_SECIS_EXTENSION
#include "SECIS/secis_neighbors.h"
#endif

/* Tons of static arrays in this one! */
static char UNUSED rcsid[] =
"$Id: barriers.c,v 1.38 2008/01/10 14:40:01 ivo Exp $";

static char *form;         /* array for configuration */
static loc_min *lmin;      /* array for local minima */

static double **rate;      /* rate matrix between basins */
static double  *dr;        /* increments to rate matrix  */

static int n_lmin;
static unsigned int max_lmin;
static int n_saddle;
static double minh;
static double energy;
    /* energy of last read structure (for check_neighbors) */
static int *POV;    /* list of last read POSET values */
static int POV_size;
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
static int print_labels = 0;
static int IS_arbitrary = 0;

static int maxlabellength = 0;

/* private functions */
static void walk_limb (hash_entry *hp, int LM, int inc, const char *tag);
static void backtrack_path_rec (int l1, int l2, const char *tag);
static int *make_sorted_index(int *truemin);
static void Sorry(char *GRAPH);
static void print_hash_entry(hash_entry *h);
static int  read_data(barrier_options opt, double *energy,char *strucb,
		      int len, int *POV);

static void merge_components(int c1, int c2);
static int comp_comps(const void *A, const void *B);

/* public functiones */
extern int cut_point;
extern char *costring(char *string);
extern char *tokenize(char *line);

int      *make_truemin(loc_min *Lmin);
loc_min  *barriers(barrier_options opt);

static int  compare(const void *a, const void *b);
void check_neighbors(void);
static void merge_basins(void);
void print_results(loc_min *L, int *tm, char *farbe);
void ps_tree(loc_min *Lmin, int *truemin, int rates);
void print_rates(int n, char *fname);
char *strip(char *s);

struct comp {
  Set *basins; /* set of basins connected by these saddles */
  char *saddle; /* one representative (first found) */
  int size;
};

static int *truecomp;
static struct comp *comp;
static int max_comp=1024, n_comp;
static int do_rates=0;
static int do_microrates=0;
static double noLP_rate=1.;
static int ligand = 0;

#define HASHSIZE (((unsigned long) 1<<HASHBITS)-1)
static hash_entry *hpool;



/* ----------------------------------------------------------- */

void set_barrier_options(barrier_options opt) {
  print_saddles = opt.print_saddles;
  bsize = opt.bsize;
  shut_up = opt.want_quiet;
  max_print = opt.max_print;
  minh = opt.minh;
  verbose = opt.want_verbose;
  print_labels = opt.label;
  switch(opt.GRAPH[0]) {
  case 'R' :    /* RNA secondary Structures */
    if (strncmp(opt.GRAPH, "RNA", 3)==0) {
      int nolp=0, shift=0, i=0;
      IS_RNA=1;
      if (opt.kT<=-300) opt.kT=37;
      kT = 0.00198717*(273.15+opt.kT);   /* kT at 37C in kcal/mol */
      move_it = RNA_move_it;
      free_move_it = RNA_free_rl;
      pack_my_structure = pack_structure;
      unpack_my_structure = unpack_structure;
      if (strstr(opt.GRAPH, "noLP")) {
	nolp=1;
	noLP_rate = opt.noLP_rate;
      }
      if(strlen(opt.MOVESET) > 0)
	switch(opt.MOVESET[0]){
	case 'S': /* Shift moves */
	  if(strncmp(opt.MOVESET,"Shift",5)==0) 
	    shift = 1;
	  break;
	case 's': /* Shift moves */
	  if(strncmp(opt.MOVESET,"shift",5)==0)
	    shift = 1;
	  break;
	case 'l': /* ligand */
	  if(strncmp(opt.MOVESET,"ligand",6)==0){
	    ligand = 1;
	    /* move_it = RNA_move_it; */
	  }
	  break;
	default:
	  fprintf(stderr, "Unknown moveset %s\n", opt.MOVESET);
	  exit (EXIT_FAILURE);
	}
      for (i=0; i < (int)strlen(opt.seq); i++){
	if (opt.seq[i] == 'T')
	  opt.seq[i] = 'U';
      }
      RNA_init(opt.seq, shift, nolp);
      if (verbose)
	fprintf(stderr, "Graph is RNA with noLP=%d, Shift=%d, ligand=%d\n", nolp, shift,ligand);
    } else Sorry(opt.GRAPH);
    break;
  case 'Q' :    /* Haming graphs */
    if (strcmp(opt.GRAPH,"Q2")==0) {   /* binary +- alphabet */
      if(strcmp(opt.MOVESET,"c")==0) {
	move_it = SPIN_complement_move_it;
	if (verbose)
	  fprintf(stderr, "Graph is Q2 with complementation moves\n");
      }
      else {
	move_it = SPIN_move_it;
	if (verbose) fprintf(stderr, "Graph is Q2\n");
      }
      pack_my_structure = pack_spin;
      unpack_my_structure = unpack_spin;
    }
    else {
      int alphabetsize=0;
      int numconv, i;
      char *ALPHA;
      ALPHA = (char *) space(strlen(opt.GRAPH)*sizeof(char));
      numconv = sscanf(opt.GRAPH,"Q%d,%s",&alphabetsize,ALPHA);
      switch(numconv) {
      case 2 :
	if((int)strlen(ALPHA)!=alphabetsize) Sorry(opt.GRAPH);
	break;
      case 1 :
	if((alphabetsize<=0)||(alphabetsize>26)) Sorry(opt.GRAPH);
	free(ALPHA);
	ALPHA = (char *) space(sizeof(char)*(alphabetsize+1));
	for(i=0;i<alphabetsize;i++) ALPHA[i] = (char) 65+i;
	break;
      default:
	Sorry(opt.GRAPH);
	break;
      }
      String_set_alpha(ALPHA);
      if(strcmp(opt.MOVESET,"c")==0){
	initialize_crankshaft();
	move_it = String_move_it_crankshaft;
	if(verbose)
	  fprintf(stderr, "Graph is Q%d with Alphabet '%s' with crankshaft moves\n",
		  alphabetsize,ALPHA);
      }
      else{
	move_it = String_move_it;
	if(verbose)
	  fprintf(stderr, "Graph is Q%d with Alphabet '%s'\n",
		  alphabetsize,ALPHA);
      }
      if(alphabetsize < 7){
	ini_pack_em(opt);
	pack_my_structure = pack_em;
	unpack_my_structure = unpack_em;
      }
      else {
	pack_my_structure = strdup;
	unpack_my_structure = strdup;
      }
    }
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
  case 'S':     /* multi objective SECIS design */
#if HAVE_SECIS_EXTENSION
    move_it = SECIS_move_it;
    pack_my_structure = strdup;
    unpack_my_structure = strdup;
#else
    fprintf(stderr,
	    "You need to reconfigure barriers with the --with-secis"
	    " option\nto use barriers SECIS design extension\n");
    exit(EXIT_FAILURE);
#endif
    break;
  case 'T' :    /* Phylogenetic Trees */
    move_it = NNI_move_it;
    pack_my_structure = strdup;
    unpack_my_structure = strdup;
    if (verbose)
      fprintf(stderr, "Graph is Trees with NNI moves\n");
    break;
  case 'X' : /* Johnson graph J(n,n/2) = balanced +/- with exchange moves */
    move_it = EXCH_move_it;
    pack_my_structure = pack_spin;
    unpack_my_structure = unpack_spin;
    break;
  case '?' : /* General graph; adjacency list on file */
    move_it = LIST_move_it;
    pack_my_structure = strdup;
    unpack_my_structure = strdup;
    IS_arbitrary = 1;
    break;
  default :
    Sorry(opt.GRAPH);
  }
  if (kT<0) {
    if (opt.kT<=-300) kT=1;
    else kT=opt.kT;
  }
  do_rates = opt.rates;
  if(opt.microrates){
    do_microrates = opt.microrates;
    do_rates = opt.microrates;
  }
}

static void Sorry(char *GRAPH) {
  fprintf(stderr,"Graph \"%s\" is not implemented\n",GRAPH);
  exit(-2);
}

static FILE *mergefile=NULL;
static int readl=0;
loc_min *barriers(barrier_options opt) {
  int length;
  double new_en=0;

  hpool = (hash_entry *) space((HASHSIZE+1)*sizeof(hash_entry));
  set_barrier_options(opt);

  length = (int) strlen(opt.seq);
  max_lmin = 16383;
  lmin = (loc_min *) space((max_lmin + 1) * sizeof(loc_min));
  n_lmin = 0;

  form = (char *) space((length+2)*sizeof(char));
  comp = (struct comp *) space((max_comp+1) * sizeof(struct comp));
  truecomp = (int *) space((max_comp+1) * sizeof(int));
  if(opt.poset) {
    POV_size = opt.poset;
    POV  = (int *) space(sizeof(int)*opt.poset);
  }
  else POV = NULL;
  ini_stapel(length);
  if (opt.ssize) {
    mergefile = fopen("saddles.txt", "w");
    if (!mergefile) fprintf(stderr, "can't open saddle file\n");
  }

  while (read_data(opt, &new_en,form,length,POV)) {
    /* fprintf(stderr, "R%s\n", form); */
    if (readl==0) mfe=energy=new_en;
    if (new_en<energy)
      nrerror("unsorted list!\n");
    if (new_en>energy) {
      /* new energy band started */
      merge_basins();
      /* fprintf(stderr, "%d %d\n", readl, lmin[1].my_pool); */
      n_comp=0;
    }
    energy = new_en;
    readl++;
    move_it(form);       /* generate all neighbor of configuration */
    /* fprintf(stderr, "M%s\n", form); */
    check_neighbors();   /* flood the energy landscape */
    reset_stapel();
    if ((n_saddle+1 == max_print) && (!opt.rates))
      break;  /* we've found all we want to know */
  }
  switch(opt.GRAPH[0]) {
  case 'Q':
    if (strcmp(opt.MOVESET,"c")==0)
      Q_mem_cleanup();
    break;
  default:
    break;
  }
  merge_basins();
  if (mergefile) fclose(mergefile);
  if(!shut_up) fprintf(stderr,
		       "read %d structures, to find %d saddles\n",
		       readl, n_saddle);

  if (max_print == 0 || max_print > n_lmin)
    max_print = n_lmin;

  lmin[0].fathers_pool = n_lmin;   /* store size here; pfs 03 2001 */
  lmin[0].E_saddle = energy + 0.001;
  lmin[0].energy = lmin[1].energy;

  if (!do_rates) {
    if (free_move_it)
      free_move_it();
    free_stapel();
  }
  free(form);
  fflush(stdout);
  if(!shut_up) fprintf(stderr, "%lu hash table collisions\n", collisions);
  free(truecomp);
  free(comp);
  return lmin;
}

int *make_truemin(loc_min *Lmin) {
  int *truemin, nlmin, i,ii;
  nlmin = Lmin[0].fathers_pool;
  truemin = (int *) space((nlmin+1)*sizeof(int));
  /* truemin[0] = nlmin; */

  for (ii=i=1; (i<=max_print)&&(ii<=n_lmin); ii++) {
    int f;
    f = lmin[ii].father;
    if (!f) lmin[ii].E_saddle = energy + 0.000001;
    if (lmin[ii].E_saddle - lmin[ii].energy - (float)minh + FLT_EPSILON >= 0.)
      truemin[ii]=i++;
    else { /* ii is not a truemin */
      lmin[f].Z += lmin[ii].Z;
      lmin[f].Zg += lmin[ii].Zg;
    }
  }
  truemin[0] = i-1;
  return truemin;
}

/*=============================================================*/

static int read_data(barrier_options opt, double *energy, char *strucb,
		     int len, int *POV){
  int   r,l;
  char *line;
  char *token;
#ifdef _DEBUG_POSET_
  static count = 1;
#endif

  line = get_line(opt.INFILE);

  if(line==NULL) return 0;
  if(strlen(line)==0) return 0;

  token=strtok(line," \t");
  if(token==NULL) return 0;
  l = strlen(token);
  if(l<1) return 0;
  if (cut_point == -1)
    strcpy(strucb,token);
  else {
    char *mystr;
    mystr = (char *) space((l+1)*sizeof(char));
    strcpy(mystr,token);
    mystr = tokenize(mystr);
    strcpy(strucb,mystr);
    free(mystr);
  }

  token = strtok(NULL," \t");
  if(token==NULL) { fprintf(stderr, "Error in input file\n"); exit(123); }
  r=sscanf(token,"%lf",energy);
  if(r<1) { fprintf(stderr, "Error in input file\n"); exit(124); }

  if(IS_arbitrary) {
    /* record the maximal length of token name for output formatting */
    if(l>maxlabellength) maxlabellength=l;
    if(l>len) {
      fprintf(stderr,"read_data():\n%s\n label too long !!\n", strucb);
      exit (111);
    }
  }
#if 0
  /*
   * removed because in the lattice protein case, the sequence is one
   * character longer than the actual SAWs
   */
  else {
    if(l != len) {
      fprintf(stderr,"read_data():\n%s\n unequal length !!\n", strucb);
      exit (112);
    }
  }
#endif

  if(opt.poset) {
    int i,x;
    for(i=0;i<opt.poset;i++) {
      token = strtok(NULL," \t");
      if(token==NULL) { fprintf(stderr, "Error in input file\n"); exit(125); }
      r = sscanf(token,"%d",&x);
      if(r!=1) { fprintf(stderr, "Error in input file\n"); exit(126); }
      POV[i]=x;
    }
#ifdef _DEBUG_POSET_
    {
      int i;
      fprintf(stderr,"POV[%4d] = {", count);
      for(i=0;i<opt.poset;i++) {
	fprintf(stderr,"%2d", POV[i]);
	if (i<opt.poset-1) fprintf(stderr,",");
      }
      fprintf(stderr, "}\n");
    }
    count++;
#endif
  }

  if(IS_arbitrary) {
    token = strtok(NULL," \t");
    if(token==NULL) put_ADJLIST(":");
    else put_ADJLIST(token);
  }

  /* the rest of the line is junk an can safely (?) be ignored */

  free(line);

  return(1);

}

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
  int basin, obasin=-1;
  hash_entry *hp, h, *down=NULL;
  Set *basins; basinT b;

  double minenergia =  100000000.0;  /* energy of lowest neighbor */
  double Zi;
  int   min_n = 1000000000;         /* index of lowest neighbor  */
  int   gradmin=0;          /* for Gradient Basins */
  int is_min=1;
  int ccomp=0;              /* which connected component */
  basins = new_set(10);

  Zi = exp((mfe-energy)/kT);

  /* foreach neighbor structure of configuration "Structure" */
  while ((p = pop())) {
    if (IS_RNA)
      if (p[strlen(p)-1] == 'D')
	p[strlen(p)-1] = '\0';
	  
    pp = pack_my_structure(p);
    h.structure = pp;


    hp = lookup_hash(&h);

    if (hp && POV_size) { /* need to check if h is dominated by hp */
      int i;
      for(i=0;i<POV_size;i++) {
	/* printf(" %d",hp->POV[i]); */
	if (POV[i] < hp->POV[i]) { hp=NULL; break; }
      }
    }
    /* check whether we've seen the structure before */
    if (hp) {
      /* because we've seen this structure before, it already */
      /* belongs to the basin of attraction of a local minimum */
      basin = hp->basin;
      if ( hp->energy < energy) is_min=0;  /* should we use hp->n here? */
      if ( hp->n < min_n ) {         /* find lowest energy neighbor */
	minenergia = hp->energy;
	min_n = hp->n;
	gradmin = hp->GradientBasin;
	down = hp;
      }

      /* careful: if input has higher precision than FLT_EPSILON
	 bad things will happen */
      if ( fabs(hp->energy - energy)<=FLT_EPSILON*fabs(energy)) {
	int tc; tc = hp->ccomp;
	while (tc != truecomp[tc]) tc = truecomp[tc];
	if (ccomp==0)
	  ccomp = tc;
	else {
	  while (ccomp != truecomp[ccomp]) ccomp = truecomp[ccomp];
	  if (ccomp != tc) merge_components(tc, ccomp);
	  ccomp = truecomp[ccomp];
	}
      }
      /* the basin of attraction of this local minimum may have been */
      /* merged with the basin of attraction of an energetically */
      /* "deeper" local minimum in a previous step */
      /* go and find this "deeper" local minimum! */
      while (lmin[basin].father) basin=lmin[basin].father;

      /* put the "deepest" local minimum into the basins-list */
      if (basin != obasin) {
	b.hp = hp; b.basin = basin;
	set_add(basins, &b);
      }
      obasin = basin;
    }
    free(pp);
  }

  /* pack read structure from subopt for putting into the hash */
  /* fprintf(stderr,"F%s\n",form); */
  pform = pack_my_structure(form);

  if (ccomp==0) {
    /* new compnent */
    Set *set;
    set = new_set(10);
    if (++n_comp>max_comp) {
      max_comp *= 2;
      comp = (struct comp*) xrealloc(comp, (max_comp+1)*sizeof(struct comp));
      truecomp = (int*) xrealloc(truecomp, (max_comp+1)*sizeof(int));
    }
    comp[n_comp].basins = set;
    comp[n_comp].saddle = pform;
    comp[n_comp].size = 0;
    truecomp[n_comp] = ccomp = n_comp;
  }

  if (is_min) {
    basinT b;
    /* Structure is a "new" local minimum */
    gradmin = ++n_lmin;        /* for Gradient Basins */
    down = NULL;
    /* need to allocate more space for the lmin-list */
    if (n_lmin > max_lmin) {
      fprintf(stderr, "increasing lmin array to %d\n",max_lmin*2);
      lmin = (loc_min *) xrealloc(lmin, (max_lmin*2+1)*sizeof(loc_min));
      memset(lmin + max_lmin +1, 0, max_lmin);
      max_lmin *= 2;
    }

    /* store configuration "Structure" in lmin-list */
    lmin[n_lmin].father = 0;
    lmin[n_lmin].structure = pform;
    lmin[n_lmin].energy = energy;
    lmin[n_lmin].my_GradPool = 0;
    lmin[n_lmin].my_pool = 1;
    lmin[n_lmin].Z = Zi;
    lmin[n_lmin].Zg = 0;
    b.basin = n_lmin; b.hp=NULL;
    set_add(basins, &b);
  }
  else comp[ccomp].size++;
  set_merge(comp[ccomp].basins, basins);

  {
    int i_lmin;
    i_lmin = (is_min) ? n_lmin : basins->data[0].basin;
    set_kill(basins);
    /* store configuration "Structure" in hash table */
    hp = hpool+readl-1;  /* (hash_entry *) space(sizeof(hash_entry)); */
    if (POV_size) {
      int i;
      hp->POV = (int *) space(sizeof(int)*POV_size);
      for(i=0;i<POV_size;i++) hp->POV[i]=POV[i];
    }

    hp->structure = pform;
    hp->energy = energy;
    hp->basin = i_lmin;
    hp->GradientBasin = gradmin;    /* for Gradient Basins */
    hp->down = down;
    hp->ccomp = ccomp;
    hp->n = readl;
    lmin[gradmin].my_GradPool++;
    lmin[gradmin].Zg += Zi;
  
    if (write_hash(hp)){
      
      hash_entry *foo = NULL;
      foo=(hash_entry*)lookup_hash(hp);
      fprintf(stderr,"%s\n",unpack_my_structure(foo->structure));
      fprintf(stderr,"%s\n",unpack_my_structure(hp->structure));
      nrerror("duplicate structure");
    }
  }

  if((is_min)&&(POV_size)) lmin[n_lmin].POV = hp->POV;
}

static void merge_basins() {
  int c, i, t;
  for (i=t=1; i<=n_comp; i++) {
    if (truecomp[i]==i)
      comp[t++]=comp[i];
    else set_kill(comp[i].basins);
  }
  n_comp = t-1;
  qsort(comp+1, n_comp, sizeof(struct comp), comp_comps);
  for (c=1; c<=n_comp; c++) { /* foreach connected component */
    /* merge all lmins connected by this component */
    static int false_lmin=0;
    int i, father, pool=0;
    double Z=0;
    basinT *basins;

    if (mergefile && (comp[c].basins->num_elem>1)) {
      const char format[2][16] = {"%13.5f %4d %s", "%6.2f %4d %s"};
      char *saddle;
      saddle = unpack_my_structure(comp[c].saddle);
      fprintf(mergefile, format[IS_RNA], energy, comp[c].size, saddle);
      free(saddle);
      for (i=0; i < comp[c].basins->num_elem; i++)
	fprintf(mergefile, " %2d", comp[c].basins->data[i].basin);
      fprintf(mergefile, "\n");
    }

    basins = comp[c].basins->data;
    father = basins[0].basin;
    while (lmin[father].father) father=lmin[father].father;

    for (i = 1; i < comp[c].basins->num_elem; i++) {
      int ii, l, r;
      ii = basins[i].basin;
      while (lmin[ii].father) ii=lmin[ii].father;
      if (ii!=father) {
	if (ii<father) {int tmp; tmp=ii; ii=father; father=tmp; l=0; r=i;}
	else {l=i; r=0;}
	/* going to merge ii with father  */
	if ((!max_print) || (ii<=max_print+false_lmin)) {
	  /* found the saddle for a basin we're gonna print */
	  if (energy-lmin[ii].energy>=minh) n_saddle++;
	  else false_lmin++;
	}

	lmin[ii].father = father;
	lmin[ii].saddle = comp[c].saddle;
	lmin[ii].E_saddle = energy;
	lmin[ii].left =  basins[l].hp;
	lmin[ii].right = basins[r].hp;
	if (bsize) {
	  lmin[ii].fathers_pool = lmin[father].my_pool;
	  pool += lmin[ii].my_pool;
	  Z += lmin[ii].Z;
	}
      }
    }
    if (bsize) {
      lmin[father].my_pool += pool + comp[c].size;
      lmin[father].Z += Z + comp[c].size * exp((mfe-energy)/kT);
    }
    set_kill(comp[c].basins);
  }
}

void mark_global(loc_min *Lmin)
{
  int i,j,k;
  int n_gmin;
  loc_min *G;

  G = (loc_min *) space( n_lmin*sizeof(loc_min) );

  Lmin[1].global = (char) 1;
  G[1] = Lmin[1];
  n_gmin = 1;

  for (i=1; i<=n_lmin; i++) {
    Lmin[i].global = (char) 1;
    for (j=1; j<=n_gmin; j++) {
      int dom;
      dom =1;
      for(k=0;k<POV_size;k++)
	if(G[j].POV[k]>=Lmin[i].POV[k]) dom=0;
      if(dom) {
	Lmin[i].global = (char) 0;
	j=n_gmin+1;
      }
    }
    if (Lmin[i].global) {
      n_gmin++;
      G[n_gmin] = Lmin[i];
    }
  }
}

/*====================*/
void print_results(loc_min *Lmin, int *truemin, char *farbe)
{
  int i,ii,j, n;
  char *struc=NULL,*laststruc=NULL;
  char *format=NULL,*formatA=NULL,*formatB=NULL;
  bool otherformat=false;

  if (POV_size) fprintf(stderr," POV_size = %d\n",POV_size);
  if (IS_arbitrary) {
    char tfor[100];
    sprintf(tfor,"%%4d %%-%ds %%6.2f %%4d %%6.2f",maxlabellength);
    format = tfor;
  }
  else if (IS_RNA){
    formatA = "%4d %s %6.2f %4d %6.2f";
    formatB = "%4d %s  %6.2f %4d %6.2f";
    format=formatA;
  }
  else
    format = "%4d %s %13.5f %4d %13.5f";
  if(verbose) printf("Using output format string '%s'\n",format);

  n_lmin = Lmin[0].fathers_pool;

  printf("     %s\n", farbe);
  for (i = 1; i <= n_lmin; i++) {
    int f;
    if ((ii = truemin[i])==0) continue;

    struc = unpack_my_structure(Lmin[i].structure);
    /* if(ligand){ */
    /*   if(laststruc != NULL){ */
    /* 	if (strncmp(laststruc,struc,strlen(struc))==0){ */
    /* 	  1;  */
    /* 	} */
    /*   } */
  
    if (cut_point > -1)
      struc = costring(struc);
    n = strlen(struc);
    f = Lmin[i].father; if (f>0) f = truemin[f];
    if(POV_size) {
      int jj;
      printf("%4d %s ", ii, struc);
   
      if(IS_RNA) printf("%6.2f ", Lmin[i].energy);
      else printf("%13.5f ", Lmin[i].energy);
      for(jj=0;jj<POV_size;jj++) printf("%6d ",Lmin[i].POV[jj]);
      if(IS_RNA) printf("%4d %6.2f",f,
			Lmin[i].E_saddle - Lmin[i].energy);
      else printf("%4d %13.5f",f,
		  Lmin[i].E_saddle - Lmin[i].energy);

      if(Lmin[i].global) printf(" *");
      else printf(" .");
    }
    else {
      if(IS_RNA && (ligand == 1)){
	if(strstr(struc,"*") == NULL){ 
	  format=formatB;
	  otherformat=true;
	}
      }
      printf(format, ii, struc, Lmin[i].energy, f,
	     Lmin[i].E_saddle - Lmin[i].energy);
      if(otherformat){
	format=formatA;
	otherformat=false;
      }
    }
    laststruc = strip(struc);
    free(laststruc);
    free(struc);

    if (print_saddles) {
      if (Lmin[i].saddle)  {
	struc = unpack_my_structure(Lmin[i].saddle);
	printf(" %s", struc);
	free(struc);
      }
      else {
	printf(" ");
	for (j=0;j<n;j++) { printf("~"); }
      }

    }
    if (bsize)
      printf (" %12ld %8ld %10.6f %8ld %10.6f",
	      Lmin[i].my_pool, Lmin[i].fathers_pool, mfe -kT*log(lmin[i].Z),
	      Lmin[i].my_GradPool, mfe -kT*log(lmin[i].Zg));
    printf("\n");
 
  }
}

/*====================*/
/* remove additional characters from structure, such as
   '*','A','B',... */
char *strip(char *s)
{
  char *p = strdup(s);
  int l = strlen(p);
  if(p[l-1] == '*') /* add more characters here is required */
    p[l-1]='\0';
  /* fprintf(stderr, "%s in \n%s out\n",s,p); */
  return p;
}

void ps_tree(loc_min *Lmin, int *truemin, int rates)
{
  nodeT *nodes;
  int i,ii;
  int nlmin;

  nlmin = Lmin[0].fathers_pool;

  if (max_print>truemin[0]) max_print=truemin[0];

  nodes = (nodeT *) space(sizeof(nodeT)*(max_print+1));
  for (i=0,ii=1; i<max_print && ii<=nlmin; ii++){
    register int s1, f;
    double E_saddle;
    if ((s1=truemin[ii])==0) continue;
    if (s1>max_print)
      nrerror("inconsistency in ps_tree, aborting");
    E_saddle = Lmin[ii].E_saddle;
    f = Lmin[ii].father;
    if (f==0) E_saddle = Lmin[0].E_saddle; /* maximum energy */

    nodes[s1-1].father = (f==0)?-1:truemin[f]-1;
    /* was truemin[f]-1; */
    if (rates) {
      double F,Ft;
      F = mfe - kT*log(Lmin[ii].Zg);
      Ft = (f>0) ? F -kT*log(rate[truemin[ii]][truemin[f]])  : E_saddle;
      nodes[s1-1].height = F;
      nodes[s1-1].saddle_height = Ft;
    }else {
      nodes[s1-1].height = Lmin[ii].energy;
      nodes[s1-1].saddle_height = E_saddle;
    }
    if (print_labels) {
      char *L;
      char *s;
      s = unpack_my_structure(Lmin[ii].structure);
      if ((POV_size)&&(Lmin[ii].global)) {
	L = (char *) space(sizeof(char)*(3+strlen(s)));
	strcat(L,s); strcat(L," *");
	nodes[s1-1].label = L;
      }
      else
	nodes[s1-1].label = strdup(s);
      free(s);
    }
    else {
      char *L=NULL, *s=NULL;
      L = (char *) space(sizeof(char)*10);
      s = unpack_my_structure(Lmin[ii].structure);
      if(s[strlen(s)-1] == '*'){
	(void) sprintf(L,"%d",s1);
	nodes[s1-1].label = L;
      }
      if((POV_size)&&(Lmin[ii].global)) {
	(void) sprintf(L,"%d *",s1);
	nodes[s1-1].label = L;
      }
      free(s);

    }
    i++;
  }
  if (rates)
    PS_tree_plot(nodes, max_print, "treeR.ps");
  else
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
static int np, max_path;

static int path_cmp(const void *a, const void *b) {
  path_entry *A, *B; int d;
  A = (path_entry *) a;
  B = (path_entry *) b;
  if ((d=strcmp(A->key, B->key))==0) return (A->num-B->num);
  else return (d);
}

/*=======*/
path_entry *backtrack_path(int l1, int l2, loc_min *LM, int *truemin) {
  int n_lmin, i, ll1=0, ll2=0;
  char *tag;
  lmin = LM;
  n_lmin = lmin[0].fathers_pool;
  for (i=1; i<=n_lmin; i++) {
    if (truemin[i]==l1) ll1=i;
    if (truemin[i]==l2) ll2=i;
  }
  np=0;
  max_path=128;
  path = (path_entry *) space(max_path*sizeof(path_entry));
  tag = (char *) space(16);
  backtrack_path_rec(ll1, ll2, tag);  path[np].hp = NULL;
  qsort(path, np, sizeof(path_entry), path_cmp);
  free(tag);
  return(path);
}

static void backtrack_path_rec (int l1, int l2, const char *tag)
{
  hash_entry h, *l1dir, *l2dir;
  int dir=1, swap=0, child, father, maxsaddle;
  /* if left==1 left points toward l2 else toward l1 */
  if (l1>l2) {
    dir = -1;
    {int t; t=l1; l1=l2; l2=t;}
  }
  child  = l2; father = l1;

  maxsaddle = child;
  /* find saddle connecting l1 and l2 */
  while (lmin[child].father != father) {
    if (lmin[child].father == 0){
      fprintf(stderr, "ERROR in backtrack_path(): ");
      fprintf(stderr,"No saddle between lmin %d and lmin %d\n", l2, l1);
      exit (1);
    }
    child = lmin[child].father;
    if (child<father) {int t; t = child; child = father; father = t;}
    if (lmin[child].E_saddle > lmin[maxsaddle].E_saddle)
      maxsaddle = child;
    /* fprintf(stderr,"f:>%d< c:>%d< %d\n", father, child, maxsaddle); */
  }
  /* found the saddle point, maxsaddle, connecting l1 and l2 */
  h.structure = lmin[maxsaddle].saddle;
  path[np].hp = lookup_hash(&h);
  strcpy(path[np].key,tag); strcat(path[np].key, "M");
  np++;
  if (np+2>=max_path) {
    max_path *= 2;
    path = (path_entry *) xrealloc(path, max_path*sizeof(path_entry));
  }

  /* which direction from saddle to l2, l1 ? */
  for (child=l2; child>0 ; child=lmin[child].father) {
    if (child==maxsaddle) {
      l2dir = lmin[maxsaddle].left;
      l1dir = lmin[maxsaddle].right;
      break;
    }
    if (child==lmin[maxsaddle].father) {
      l2dir = lmin[maxsaddle].right;
      l1dir = lmin[maxsaddle].left;
      break;
    }
  }

  /* branch to l2,  else saddle==l2 and we're done */
  if (l2dir)
    walk_limb (l2dir, l2, -dir, tag);
  /* branch to l1 (to father) */
  if (l1dir) /* else saddle==l1 and we're done */
      walk_limb (l1dir, l1, dir, tag);
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
    if (np+2>=max_path) {
      max_path *= 2;
      path = (path_entry *) xrealloc(path, max_path*sizeof(path_entry));
    }
  }

  if (inc<0) tmp[strlen(tmp)-1] = '\0';

  /* fprintf(stderr, "walk towards %d lands in %d\n", LM, htmp->basin); */

  /* wrong local minimum start cruising again */
  if (htmp->basin != LM) {
    if (inc == -1)
      backtrack_path_rec (htmp->basin, LM, tmp);
    else
      backtrack_path_rec (LM, htmp->basin, tmp);
  }
  free(tmp);
}

void print_path(FILE *PATH, path_entry *path, int *tm) {
  int i;
  for (i=0; path[i].hp; i++) {
    char c[6] = {0,0,0,0}, *struc;
    if (path[i].hp->down==NULL) {
      sprintf(c, "L%04d", tm[path[i].hp->basin]);
    } else
      if (path[i].key[strlen(path[i].key)-1] == 'M')
	c[0] = 'S';
      else c[0] = 'I';
    struc = unpack_my_structure(path[i].hp->structure);
    fprintf(PATH, "%s (%6.2f) %-5s\n", struc,  path[i].hp->energy, c);
    free(struc);
  }
}

static void merge_components(int c1, int c2) {
  if (comp[c1].size<comp[c2].size) {int cc; cc=c1; c1=c2; c2=cc;}
  comp[c1].size += comp[c2].size;
  truecomp[c2]=c1;
  set_merge(comp[c1].basins, comp[c2].basins);
}

static int comp_comps(const void *A, const void *B) {
  struct comp *a, *b;
  int r, i=0;
  a = (struct comp *)A;
  b = (struct comp *)B;
  for (i=0; i<a->basins->num_elem && i<b->basins->num_elem; i++) {
    r = a->basins->data[i].basin - b->basins->data[i].basin;
    if (r!=0) return r;
  }
  return (i==a->basins->num_elem)? -1:1;
}

static void print_hash_entry(hash_entry *h) {
  int down=0;
  if (h->down) down=h->down->n;
  fprintf(stderr, "%2d %s %6.2f %2d %2d %2d %2d\n", h->n, h->structure,
	 h->energy, h->basin, h->GradientBasin, h->ccomp, down);
}

map_struc get_mapstruc(char *p, loc_min *LM, int *tm) {
  hash_entry *hp, h;
  char *pp, *struc;
  int min, gradmin, tmin, tgradmin;
  map_struc ms;
  
  pp = pack_my_structure(p);
  h.structure = pp;
  hp = lookup_hash(&h);

  if (hp==NULL) {
    fprintf(stderr, "get_mapstruc: structure not in hash\n");
    return ms;
  }
	    
  min = hp->basin;
  while (tm[min]==0) {
    min = LM[min].father;
  }
  gradmin = hp->GradientBasin;
  while (tm[gradmin]==0) {
    gradmin = LM[gradmin].father;
  }

  if (gradmin == 0) {
    fprintf(stderr, "get_mapstruc: gradient minimum not yet assigned\n");
    return ms;
  }

  ms.structure = unpack_my_structure(LM[gradmin].structure);
  ms.n = hp->n;
  ms.energy = hp->energy;
  ms.min = min;
  ms.truemin = tm[min];
  ms.gradmin = gradmin;
  ms.truegradmin =  tm[gradmin];
  free(pp);
  return (ms);
}


void print_rates(int n, char *fname) {
  int i,j;
  FILE *OUT;
#define BINRATES
#ifdef BINRATES
  FILE *BINOUT;
  char *binfile = "rates.bin";
  double tmprate;
  BINOUT = fopen(binfile, "w");
  if (!BINOUT){
    fprintf(stderr, "could not open file pointer 4 binary outfile\n");
    exit(101);
  }
  /* first write dim to file */
  fwrite(&n,sizeof(int),1,BINOUT);
  for(i=1;i<=n;i++)
    for(j=1;j<=n;j++){
      tmprate = rate[j][i];
      fwrite(&tmprate,sizeof(double),1,BINOUT);
    }
  fprintf(stderr, "rate matrix written to binfile\n");
  fclose(BINOUT);
#endif

  OUT = fopen(fname, "w");
  if (!OUT) {
    fprintf(stderr, "could not open rates file %s for output\n", fname);
    return;
  }

  for (i=1; i<=n; i++) {
    for (j=1; j<=n; j++)
      fprintf(OUT, "%10.4g ", rate[i][j]);
    fprintf(OUT, "\n");
    free(rate[i]);
  }
  free(rate);
  fclose(OUT);
}

void compute_rates(int *truemin, char *farbe) {
  int i, j, ii, r, gb, gradmin,n, rc, *realnr;
  char *p, *pp, *form, newsub[10]="new.sub", mr[15]="microrates.out";
  hash_entry *hpr, h, *hp;
  double Zi;
  FILE *NEWSUB=NULL, *MR=NULL;;

  if(ligand==1)
    move_it = RNA_move_it_rates;
  n = truemin[0];
  rate = (double **) space((n + 1) * sizeof(double *));
  dr   = (double  *) space((n + 1) * sizeof(double));
  for (i=1; i<=n; i++)
    rate[i] = (double *) space((n + 1) * sizeof(double));
  if(do_microrates){
    realnr = (int *)space((readl+1) * sizeof(int));
    MR = fopen(mr, "w");
    NEWSUB = fopen(newsub, "w");
    fprintf(NEWSUB, "%s %6.2f\n", farbe, 100*mfe);
    fflush(NEWSUB);
    fprintf(MR, ">%d states\n", readl);
  }

  for (rc=1, r=0; r<readl; r++) {
    int b;
    hpr= &hpool[r];
    Zi = exp((mfe-hpr->energy)/kT);
    gradmin = hpr->GradientBasin;
    while (truemin[gradmin]==0) gradmin = lmin[gradmin].father;
    gradmin=truemin[gradmin];
    if (gradmin>n) continue;
    for (b=hpr->basin; b>1; b=lmin[b].father);
    form = unpack_my_structure(hpr->structure);
    move_it(form);       /* generate all neighbors of configuration */

    for (i=0; i<=n; i++) dr[i]=0;
    while ((p = pop())) {
      int double_move = 0;
      if (IS_RNA) {
	if (p[strlen(p)-1]=='D') {
	  double_move = 1;
	  p[strlen(p)-1]='\0';
	}
	else
	  double_move=0;
      }
	  
      pp = pack_my_structure(p);
      h.structure = pp;
      /* check whether we've seen the structure before */
      if ((hp = lookup_hash(&h)))
	if (hp->n <=r ) {
	  gb = hp->GradientBasin;
	  while (truemin[gb]==0) gb = lmin[gb].father;
	  gb = truemin[gb];
	  if (gb<=n)
	    dr[gb] += (double_move)?(noLP_rate*Zi):Zi;
	  if (do_microrates && b) {
	    double rate,dg;
	    dg = hpr->energy - hp->energy;
	    rate = exp(-dg/kT);
	    fprintf(MR,"%10d %8d %15.12f 1\n",rc,realnr[hp->n],rate);
	  }
	}
      free(pp);
    }
    if (do_microrates && b){
      fprintf(NEWSUB, "%s %6.2f %i %i\n", form, hpr->energy, gradmin, hpr->basin);
      fflush(NEWSUB);
      realnr[hpr->n]=rc++;
    }
    for (i=1; i<=n; i++) {
      rate[i][gradmin] += dr[i];
      rate[gradmin][i] += dr[i];
    }
    free(form);
    reset_stapel();
  }

  fprintf(stderr, "done with 2nd pass\n" );
  free(dr);

  for (i=ii=1; i<=n; i++, ii++) {
    while (truemin[ii]!=i) ii++;
    for (j=1; j<=n; j++)
      rate[i][j] /= lmin[ii].Zg;
  }
  if(do_microrates){
    free(realnr);
    fclose(NEWSUB);
    fclose(MR);
  }
  if (free_move_it)
    free_move_it();
  free_stapel();
}
