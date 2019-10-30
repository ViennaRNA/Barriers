/* trees.h */

#ifndef BARRIERS_TREES_H
#define BARRIERS_TREES_H

#include "tree_types.h"

Tree *
PTAlloc(int n);

Tree *
Make3Tree(void);

int
Numtree(int N);

Tree *
MakeAllNplus1Trees(Tree *TT);

Tree *
MakeRandomTree(int n);

Tree **
MakeTreesUpTo(int N);

void
PrintTree(Tree T);

void
FreeTree(Tree *T,
         int  n);

void
FreeAllTreesUpTo(int  N,
                 Tree **T);

Tree *
NNI_Move(Tree TT,
         int  edge,
         int  which);

Tree *
Make_all_NNI(Tree TT);

Tree *
MakeRandomNNI(Tree TT);

char *
Tree2string(Tree t);

char *
subtrees(Tree T,
         int  edge,
         int  vertex,
         int  *smallest_leaf);

int
number_of_leaves(char *st);

Tree *
string2Tree(char *s);

int
fill_T(Tree *T,
       char *s,
       int  call);

char *
InteriorTreeString(char *s);

#endif

/* End of file */
