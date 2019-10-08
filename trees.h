/* trees.h */
#include "tree_types.h"

extern int Numtree(int N);


extern Tree *PTAlloc(int n);


extern Tree *Make3Tree(void);


extern Tree *MakeAllNplus1Trees(Tree *TT);


extern Tree *MakeRandomTree(int n);


extern Tree **MakeTreesUpTo(int N);


extern Tree *NNI_Move(Tree  TT,
                      int   edge,
                      int   which);


extern Tree *Make_all_NNI(Tree TT);


extern Tree *MakeRandomNNI(Tree TT);


extern void PrintTree(Tree t);


extern void FreeTree(Tree *T,
                     int  n);


extern void FreeAllTreesUpTo(int  N,
                             Tree **T);


extern char *Tree2string(Tree t);


extern Tree *string2Tree(char *s);


extern int  number_of_leaves(char *st);


extern char *InteriorTreeString(char *s);


/* End of file */
