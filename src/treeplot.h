#ifndef BARRIERS_TREEPLOT_H
#define BARRIERS_TREEPLOT_H

/* treeplot.h */
typedef struct node {
  float height;         /* height (energy, time, whatever) of this leaf    */
  float saddle_height;  /* height of internal node that connects this leaf */
  unsigned long father;         /* node with which it connects                     */
  char  *label;         /* label string, if NULL use index+1               */
} nodeT;

void PS_tree_plot(nodeT *nodes,
                  unsigned long  n,
                  char  *filename);


#endif
/* End of file */
