/* treeplot.h */
typedef struct node {
  float height;         /* height (energy, time, whatever) of this leaf    */
  float saddle_height;  /* height of internal node that connects this leaf */
  int father;           /* node with which it connects                     */
  char *label;          /* label string, if NULL use index+1               */ 
} nodeT;

void PS_tree_plot(nodeT *nodes, int n, char *filename);

/* End of file */
