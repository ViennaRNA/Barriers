/* Last changed Time-stamp: <2001-04-08 14:49:28 ivo> */
/* treeplot.c */
/* modified version from ViennaRNA-package */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "system.h"
#include "utils.h"

static char UNUSED rcsid[]= "$Id: treeplot.c,v 1.2 2001/04/09 08:02:00 ivo Exp $";

typedef struct node {
  float height;         /* height (energy, time, whatever) of this leaf   */
  float saddle_height;  /* height of internal node that connects this leaf */
  int father;           /* node with which it connects                    */
  char *label;          /* label string, if NULL use number+1             */ 
} nodeT;

static nodeT *leafs;
static int cmp_saddle(const void *, const void*);

typedef struct link {
  int x;
  struct link *next;
} linkT;

void PS_tree_plot(nodeT *nodes, int n, char *filename) {
  /* plot tree with leaf-nodes nodes, to file filename (stdout if NULL) */
  FILE *out;
  int i, k, f=0, *sindex;
  linkT *chain, *l;
  int bbox[4] = {72, 144, 522, 700};
  if (filename) {
    out = fopen(filename, "w");
    if (!out) {
      fprintf(stderr, "can't open file %s, aborting plot\n", filename);
      return;
    }
  } else
    out = stdout;

  /* make index list, sorted by saddle height */
  sindex = (int *) space(n * sizeof(int));
  for (i=0; i<n; i++) sindex[i]=i;
  leafs =  nodes;
  qsort(sindex, n, sizeof(int), cmp_saddle);

  /* make order (x-coordinates) of leafs, we use a linked list for this */
  chain = (linkT *) space(n * sizeof(linkT));
  /* start connecting the links of the chain */
  for (i=0; i<n; i++) { /* n-1 merges */
    k = sindex[i]; f = nodes[k].father; /* ith saddle merges k and f */
    if (k==f) continue;  /* lowest node doesn't merge */
    for (l= &chain[f]; l->next!=NULL; l = l->next); 
    l->next=&chain[k]; /* attach child to chain of father */
  }
  /* chain[f] now starts the ordered chain, next fill in the num field */
  for (i=0, l= &chain[f]; l!=NULL; l = l->next) l->x = i++;

  fprintf(out, "%%!PS-Adobe-2.0 EPSF-1.2\n"
	  "%%%%Title: TreePlot\n"
	  "%%%%Creator: treeplot.c\n"
	  "%%%%CreationDate: %s", time_stamp());
  fprintf(out, "%%%%Pages: 1\n"
	  "%%%%BoundingBox: %d %d %d %d\n",
	  bbox[0], bbox[1], bbox[2], bbox[3]);
  fprintf(out, "/treedict 100 dict def\n"
	  "treedict begin\n"
	  "  /cmtx matrix currentmatrix def\n"
	  "  /STR 128 string def\n"
	  "%% - => -\n"
	  "  /Labeloffset {\n"
	  "     /Lo [\n"
	  "      (X) stringwidth pop %% width\n"
	  "      newpath 0 0 moveto\n"
	  "      (X) true charpath\n"
	  "      flattenpath pathbbox\n"
	  "      pop exch pop exch sub neg 2 div %% height\n"
	  "     ] def\n"
	  "  } def\n"
	  "%% str => -\n"
	  "  /Rotshow {\n"
	  "    gsave\n"
	  "     cmtx setmatrix 90 rotate\n"
	  "     dup stringwidth pop \n"
	  "     Lo aload pop\n"
	  "     3 1 roll add neg exch\n"
	  "      rmoveto show\n"
	  "    grestore\n"
	  "  } def\n"
	  "%% - => -\n"
	  "  /Drawlabels {\n"
	  "   0 LEAF {\n"
	  "      aload pop moveto\n"
	  "      dup LABEL exch get STR cvs Rotshow\n"
	  "      1 add\n"
	  "    } forall pop\n"
	  "  } def\n"
	  "%% - => -\n"
	  "  /Connectlmins {\n"
	  "    newpath\n"
	  "    SADDEL {\n"
	  "      aload pop  %% => child father height\n"
	  "      LEAF 4 -1 roll get aload pop %% => f h cx cy\n"
	  "      1 index exch moveto          %% => f h cx\n"
	  "      2 copy exch lineto           %% => f h cx\n"
	  "      LEAF 3 index get aload pop   %% => f h cx fx fy\n"
	  "      1 index dup 5 index lineto   %% => f h cx fx fy fx\n"
	  "      exch lineto\n"
	  "      add 2 div exch [ 3 1 roll]   %% => father [nx h]\n"
	  "      LEAF 3 1 roll put\n"
	  "    } forall\n"
	  "    gsave\n"
	  "      cmtx setmatrix stroke\n"
	  "    grestore\n"
	  "  } def\n"
	  "%% data starts here!!!\n"
	  "  /LABEL [");

  /* print label array */
  for (i=0; i<n; i++) {
    if (i%10 == 0)  fprintf(out, "\n   ");
    if (nodes[i].label) fprintf(out, "(%s) ", nodes[i].label);
    else fprintf(out, "%3d ", i+1);
  }
  fprintf(out, "\n  ] def\n");

  /* print leaf node coordinates */
  fprintf(out, "%% leaf node coordinates\n"
	  "  /LEAF [");
  for (i=0; i<n; i++) {
    if (i%5 == 0)  fprintf(out, "\n   ");
    fprintf(out, "[%-3d %7.3f] ", chain[i].x, nodes[i].height);
  }
  fprintf(out, "  \n] def\n");

  /* print internal node coordinates */
  fprintf(out, "%% internal nodes (saddle) coordinates, sorted by height\n"
	  "  /SADDEL [");
  for (i=0; i<n-1; i++) {
    k=sindex[i];
    if (i%4 == 0)  fprintf(out, "\n   ");
    fprintf(out, "[%3d %3d %7.3f] ",k,nodes[k].father, nodes[k].saddle_height);
  }
  free(chain);  free(sindex);
  fprintf(out, "  \n] def\n"
	       "end\n");
  fprintf(out, "%%%%EndProlog\n"
	  "treedict begin\n"
	  "  /fsize 8 def\n"
	  "  /Helvetica findfont\n"
	  "  fsize scalefont setfont\n"
	  "  Labeloffset\n"
	  "  %d %d fsize add translate\n", bbox[2], bbox[1]);
  fprintf(out, "  %d %d sub LEAF length div %% x-scale\n", bbox[0], bbox[2]);
  fprintf(out, "  %d %d fsize add sub\n", bbox[3]-1, bbox[1]);
  fprintf(out, "  SADDEL dup length 1 sub get 2 get %% max_height\n"
	  "LEAF 0 get 1 get sub %% energy interval\n"
	  "  div scale\n"
	  "  .5 LEAF 0 get 1 get neg translate\n"
	  "  Drawlabels\n"
	  "  Connectlmins\n"
	  "  showpage\n"
	  "end\n"
	  "%%%%EOF\n");
  if (filename) fclose(out);
}

static int cmp_saddle(const void *A, const void *B) {
  float diff;
  diff = leafs[*((int *)A)].saddle_height - 
    leafs[*((int *)B)].saddle_height;
/*    fprintf(stderr, "%d %d %f\n",*((int *)A), *((int *)B), diff);  */
  if (diff < -1e-6) return -1;
  else if (diff>1e-6) return 1;
  diff = leafs[*((int *)A)].height - leafs[*((int *)B)].height;
  return (diff<0)?-1:1;
}
  
