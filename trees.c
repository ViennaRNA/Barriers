/* trees.c */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tree_types.h"
#include "utils.h"

#ifndef DEBUG
#define  DEBUG 0
#endif
void PrintTree(Tree t);
int Numtree(int N);
Tree *PTAlloc(int n);
Tree *Make3Tree(void);
Tree *MakeAllNplus1Trees(Tree *TT);
Tree *Make_all_NNI(Tree TT);
Tree *NNI_Move(Tree TT, int edge, int which);
Tree *MakeRandomNNI(Tree TT);
void FreeTree(Tree *T,int n);
Tree **MakeTreesUpTo(int N);
Tree *MakeRandomTree(int N);

char *Tree2string(Tree t);
char *subtrees(Tree T, int edge, int vertex, int *smallest_leaf);
int number_of_leaves(char *st);
Tree *string2Tree(char *s);
int fill_T(Tree *T, char *s,int call);
     
Tree *PTAlloc(int n)
{
  Tree *T;
  T = (Tree *) space(sizeof(Tree));
  T->size = n;
  T->L = (Node *) space((n+1)*sizeof(Node));
  T->I = (Node *) space((n-1)*sizeof(Node));
  T->E = (Edge *) space((2*n-2)*sizeof(Edge));
  return T;
}

Tree *Make3Tree(void)
{
  Tree *T;
  T = PTAlloc(3);
  T->L[1].e1 = 1;
  T->L[2].e1 = 2;
  T->L[3].e1 = 3;
  T->I[1].e1 = 1;
  T->I[1].e2 = 2;
  T->I[1].e3 = 3;
  T->E[1].v1 = 1; T->E[1].v2 = -1;
  T->E[2].v1 = 2; T->E[2].v2 = -1;  
  T->E[3].v1 = 3; T->E[3].v2 = -1;  
  return T;
}

int Numtree(int N)
{
  int i,k;
  if (N<1) return 0;
  if (N<=3) return 1;
  k=1;
  for (i=3;i<=N;i++) k*=(2*i-5);
  return k;
}

Tree *MakeAllNplus1Trees(Tree *TT)
{
  int i,j,k,l,N;

  Tree *T;
  Tree *TempTree;
  int tmpv;

  N = TT->size;
  T = (Tree*) space((Numtree(N+1))*sizeof(Tree));

  for (i=0;i<Numtree(N);i++) { 
    for (j=1;j<=2*N-3;j++) {
      k = (2*N-3)*i+(j-1);

      /* copy tree i */
      TempTree = PTAlloc(N+1); T[k] = TempTree[0];

      for(l=1;l<=N;  l++) T[k].L[l].e1=TT[i].L[l].e1;
	
      for(l=1;l<=N-2;l++) {
	T[k].I[l].e1 = TT[i].I[l].e1;
	T[k].I[l].e2 = TT[i].I[l].e2;	
	T[k].I[l].e3 = TT[i].I[l].e3;  
      }

      for(l=1;l<=2*N-3;l++) {
	T[k].E[l].v1 = TT[i].E[l].v1;
	T[k].E[l].v2 = TT[i].E[l].v2;
      }

      /* now replace edge j with an interior vertex and the extra
	 leave */

      tmpv =  T[k].E[j].v2;
      T[k].E[j].v2   = -(N-1);        /* new interior vertex */
      if (tmpv > 0 )  {    /* leave */
	T[k].L[tmpv].e1 = 2*N-3+1;
      }
      else {               /* interior node */
        if (T[k].I[-tmpv].e1==j) T[k].I[-tmpv].e1 = 2*N-3+1;
	else if (T[k].I[-tmpv].e2==j) T[k].I[-tmpv].e2 = 2*N-3+1;
	else if (T[k].I[-tmpv].e3==j) T[k].I[-tmpv].e3 = 2*N-3+1;
      }
      
      T[k].E[2*N-3+1].v1 = -(N-1);
      T[k].E[2*N-3+1].v2 = tmpv;
      T[k].I[N-1].e1 = j;
      T[k].I[N-1].e2 = 2*N-3+1;
      T[k].I[N-1].e3 = 2*N-3+2;
      T[k].L[N+1].e1 = 2*N-3+2;   
      T[k].E[2*N-3+2].v1 = N+1;    
      T[k].E[2*N-3+2].v2 = -(N-1);

      free(TempTree);
    }
  }
  return T;
}

Tree *MakeRandomTree(int n)
{
  int j,l,N;
  int tmpv;
  Tree *T,*TT;
  
  if (n<3) return NULL;

  T = Make3Tree();

  for (N=3; N<n; N++) {
    (void) urn(); /* what's this ??? */
    j = int_urn(1,2*N-3);  /* insert at edge j */
    
    /* copy tree */
    TT = PTAlloc(N+1); 
    
    for(l=1;l<=N;  l++) TT->L[l].e1=T->L[l].e1;
    
    for(l=1;l<=N-2;l++) {
      TT->I[l].e1 = T->I[l].e1;
      TT->I[l].e2 = T->I[l].e2;	
      TT->I[l].e3 = T->I[l].e3;  
    }
    
    for(l=1;l<=2*N-3;l++) {
      TT->E[l].v1 = T->E[l].v1;
      TT->E[l].v2 = T->E[l].v2;
    }

    FreeTree(T,1);
    
    /* now replace edge j with an interior vertex and the extra
       leave */
    
    tmpv =  TT->E[j].v2;
    TT->E[j].v2   = -(N-1);        /* new interior vertex */
    if (tmpv > 0 )  {    /* leave */
      TT->L[tmpv].e1 = 2*N-3+1;
    }
    else {               /* interior node */
      if (TT->I[-tmpv].e1==j) TT->I[-tmpv].e1 = 2*N-3+1;
      else if (TT->I[-tmpv].e2==j) TT->I[-tmpv].e2 = 2*N-3+1;
      else if (TT->I[-tmpv].e3==j) TT->I[-tmpv].e3 = 2*N-3+1;
    }
    
    TT->E[2*N-3+1].v1 = -(N-1);
    TT->E[2*N-3+1].v2 = tmpv;
    TT->I[N-1].e1 = j;
    TT->I[N-1].e2 = 2*N-3+1;
    TT->I[N-1].e3 = 2*N-3+2;
    TT->L[N+1].e1 = 2*N-3+2;   
    TT->E[2*N-3+2].v1 = N+1;    
    TT->E[2*N-3+2].v2 = -(N-1);
    
    T=TT;
  }

  return TT;
}

Tree **MakeTreesUpTo(int N)
{
  int i;
  Tree **T;

  if (N<3) return NULL;
  T = (Tree **) space((N+1)*sizeof(Tree*));
  T[0] = NULL;
  T[1] = NULL;
  T[2] = NULL;
  for (i=3; i<=N; i++) T[i] = space( Numtree(i)*sizeof(Tree *));

  T[3] = Make3Tree();
  for (i=4; i<=N; i++) T[i] = MakeAllNplus1Trees(T[i-1]);
  
  return T;
}

void PrintTree(Tree T)
{
  int i,n;
  n = T.size;
  printf("\nSize: %d\n",T.size);
  printf("Leaves: \n");
  for (i=1;i<=n; i++) printf("%d: %d\n",i,T.L[i].e1);
  printf("Interior: \n");
  for (i=1;i<=n-2;i++)
    printf("%d: %d %d %d\n",i,T.I[i].e1,T.I[i].e2,T.I[i].e3);
  printf("Edges:\n");
  for (i=1;i<=2*n-3;i++)
    printf("%d: %d %d\n",i,T.E[i].v1,T.E[i].v2);
  printf("\n");
}

void FreeTree(Tree *T,int n)
{
  int i;
  for(i=0;i<n;i++) {
    free(T[i].L);
    free(T[i].I);
    free(T[i].E);
  }
  free(T);
}

void FreeAllTreesUpTo(int N, Tree **T)
{
  int i,n;
  for (i=3;i<=N;i++) {
    n = Numtree(i);
    FreeTree(T[i],n);
  }
  free(T);
}

Tree *NNI_Move(Tree TT, int edge, int which)
{
  int v1,v2,e1,e2,etmp;
  Tree *T;
  int N,l;

  if (DEBUG) printf("In NNI Move\n");
  if (DEBUG) PrintTree(TT);

  v1 = TT.E[edge].v1;
  v2 = TT.E[edge].v2;
  if (DEBUG) printf("&&& %d: %d %d\n",edge, v1, v2);
  if ((v1>0)||(v2>0)) return NULL;   /* null pointer if not interior edge */
  
  N = TT.size;
  
  T = PTAlloc(N); 

  for(l=1;l<=N;  l++) T->L[l].e1=TT.L[l].e1;
  
  for(l=1;l<=N-2;l++) {
    T->I[l].e1 = TT.I[l].e1;
    T->I[l].e2 = TT.I[l].e2;	
    T->I[l].e3 = TT.I[l].e3;  
  }
  
  for(l=1;l<=2*N-3;l++) {
    T->E[l].v1 = TT.E[l].v1;
    T->E[l].v2 = TT.E[l].v2;
  }

  if (DEBUG) PrintTree(T[0]);
  
  /* reorder edges so that "edge" is e3 */

  if (edge == T->I[-v1].e1) {
    etmp = T->I[-v1].e3; T->I[-v1].e3 = edge; T->I[-v1].e1 = etmp; }
  else if (edge == T->I[-v1].e2) {
    etmp = T->I[-v1].e3; T->I[-v1].e3 = edge; T->I[-v1].e2 = etmp; }
  if (edge == T->I[-v2].e1) {
    etmp = T->I[-v2].e3; T->I[-v2].e3 = edge; T->I[-v2].e1 = etmp; }
  else if (edge == T->I[-v2].e2) {
    etmp = T->I[-v2].e3; T->I[-v2].e3 = edge; T->I[-v2].e2 = etmp; }

  /* identify the edges that are to be swapped */
  
  e1 = T->I[-v1].e2;
  if (which == 1)  e2 = T->I[-v2].e1;
  else             e2 = T->I[-v2].e2;

  if (DEBUG) printf("** %d %d \n",e1,e2);

  /* make the vertices on "edge" v2 on both edges */
  if (T->E[e1].v1 == v1) { T->E[e1].v1 = T->E[e1].v2; T->E[e1].v2=v1; }
  if (T->E[e1].v1 == v2) { T->E[e1].v1 = T->E[e1].v2; T->E[e1].v2=v2; }
  if (T->E[e2].v1 == v1) { T->E[e2].v1 = T->E[e2].v2; T->E[e2].v2=v1; }
  if (T->E[e2].v1 == v2) { T->E[e2].v1 = T->E[e2].v2; T->E[e2].v2=v2; }

  if (DEBUG) PrintTree(T[0]);

  /* and now swap */

  T->E[e1].v2  = v2;
  T->E[e2].v2  = v1;
  T->I[-v1].e2 = e2;
  if (which == 1)  T->I[-v2].e1 = e1;
  else             T->I[-v2].e2 = e1;
  
  if (DEBUG) PrintTree(T[0]);
 
  return T;
}


Tree *Make_all_NNI(Tree TT)
{
  int k,l,n;
  Tree *NL;
  Tree *TempTree;
  n = TT.size;
  if(DEBUG) printf("Ausgangsbaum\n");
  if(DEBUG) PrintTree(TT);
  NL=(Tree*) space(2*(n-3)*sizeof(Tree));
  if(DEBUG) printf("NNI Neighbors:\n");
  k=0;
  for (l=1; l<=2*n-3; l++) {
    if(DEBUG) printf(">> l=%d\n",l);
    TempTree = NNI_Move(TT,l,1);
    if (TempTree!=NULL) {
      NL[k] = TempTree[0]; 
      if(DEBUG) printf(">> k=%d\n",k);
      if(DEBUG) PrintTree(NL[k]);
      k++;
      free(TempTree);
      TempTree = NNI_Move(TT,l,2);
      NL[k] = TempTree[0]; 
      if(DEBUG) printf(">> k=%d\n",k);
      if(DEBUG) PrintTree(NL[k]);
      k++;
      free(TempTree);
    }
  }

  return NL;
}

Tree *MakeRandomNNI(Tree TT)
{
  int i,k,l,n;
  Tree *NL;

  n = TT.size;
  
  k = int_urn(1,n-3);
   /* find the index l of the k-th *interior* edge */              
  i = 0; l=0;
  while (i<k) {
    l++;
    if( (TT.E[l].v1<0) && (TT.E[l].v2<0) ) i++ ;
  }
  i= (urn()>0.5)?1:0;

  NL = NNI_Move(TT,l,i);
  
  return NL;
}

/* ----------------------------------------------------------------- */

char *Tree2string(Tree t)
{
  int k,n,len;
  int e1,v1;
  char *s;
  char *s1;

  n=t.size;
  
  e1 = t.L[1].e1;
  if(t.E[e1].v1==1) v1 = t.E[e1].v2;
  else v1 = t.E[e1].v1;

  s1 = subtrees(t,e1,v1,&k);
  len = strlen(s1)+5;;
  s = (char *) space(sizeof(char)*(len+1)); s[0]='\0';
  strcat(s,"((1)");
  strcat(s,s1);
  strcat(s,")");
  free(s1);
  return s;
  
}
  
char *subtrees(Tree T, int edge, int vertex, int *smallest_leaf)
{
  int to1e, to2e, to1v, to2v;
  int len;
  char *s, *s1, *s2, ss[20];
  int sl1,sl2;
  
  if(vertex>0) {
    sprintf(ss,"%d",vertex);
    len = strlen(ss);
    s = (char *) space(sizeof(char)*(len+1));
    strcpy(s,ss);
    *smallest_leaf = vertex;
    return s;
  }

  if (T.I[-vertex].e1==edge) {
    to1e=T.I[-vertex].e2; to2e=T.I[-vertex].e3;
  }
  else {
    to1e=T.I[-vertex].e1;
    if(T.I[-vertex].e2==edge) 
      to2e = T.I[-vertex].e3;
    else
      to2e = T.I[-vertex].e2;
  }
  if (T.E[to1e].v1==vertex) 
    to1v = T.E[to1e].v2;
  else to1v = T.E[to1e].v1;
  if (T.E[to2e].v1==vertex) 
    to2v = T.E[to2e].v2;
  else to2v = T.E[to2e].v1; 

  s1 = subtrees(T,to1e,to1v,&sl1);
  s2 = subtrees(T,to2e,to2v,&sl2);

  *smallest_leaf = sl1;
  if(sl2<sl1) { s=s1; s1=s2; s2=s; *smallest_leaf =sl2;}

  len = 4+strlen(s1)+strlen(s2);

  s = (char *) space(sizeof(char)*(len+1));
  strcat(s,"("); strcat(s,s1); strcat(s,")("); strcat(s,s2); strcat(s,")");

  free(s1);free(s2);
  return s;
}
	
/* --------------------------------------------------------- */

int number_of_leaves(char *st)
{
  int i,k,l;
  char last;
  l = strlen(st);

  k=0;
  last =st[0];
  for(i=1;i<l;i++) {
    if(st[i]==')') {
      if (last=='(') k++;
      last =')';
    }
    if(st[i]=='(') last = '(';
  }
  return k;
}

/* --------------------------------------------------------- */

Tree *string2Tree(char *s) {
  int i,j,k,n,len,e1,e2;
  Tree *T;
  char *s1,*s2;

  len = strlen(s);
  n = number_of_leaves(s);
  if (n==3) {
    T = Make3Tree();
    return T;
  }
  T = PTAlloc(n);
  /* printf("size %d\n",T->size); */
  if (n==1) return T;
  if (n==2) {
    T->L[1].e1 = 1;
    T->L[2].e1 = 1;
    T->E[1].v1 = 1;
    T->E[1].v2 = 1;
    return T;
  }
  /* consider the case n>4 */

  for(i=1;i<=n;i++) {
    T->L[i].e1 = i;
    T->E[i].v1 = i;
  }

  k=0;
  for(i=4;i<len-1;i++) {
    if( s[i]=='(' ) k++;
    if( s[i]==')' ) k--;
    if(k==0) break;
  }
  j=i;

  s1 = (char *) space(sizeof(char)*(i-4+1+1));
  for(i=0;i<j-3;i++) s1[i]=s[i+4];
  s1[i]='\0';

  s2 = (char *) space(sizeof(char)*(len-(i-3)));
  for(i=j+1;i<len-1;i++) s2[i-(j+1)] = s[i];
  s2[i-(j+1)]='\0';

  e1 = fill_T(T,s1,1); /*initialize at first call */
  e2 = fill_T(T,s2,0);

  free(s1); free(s2);

  T->E[1].v2 = -1;  
  T->E[e2].v2 = -1;  
  T->E[e1].v2 = -1;  
  T->I[1].e1 = e1;
  T->I[1].e2 = e2;
  T->I[1].e3 = 1;

  return T;
}

int fill_T(Tree *T, char *s, int call)
{
  static int int_vertex;
  static int int_edge;
  
  int leaf;
  int e1, e2,n;
  int this_vertex;
  int i,j,k;
  int len;
  char *s1, *s2;
  
  n = T->size;

  if(call) {    /* dirty trick for initialization */
    int_edge   = n;
    int_vertex = 1;
  }
    
  len = strlen(s);
  if( (s[0]!='(')||(s[len-1]!=')') ) nrerror("This should not happen");

  if(s[1]!='(') {   /* leaf */
    sscanf(s,"(%d)",&leaf);
    T->E[leaf].v1=leaf;
    /* printf("> %d *\n",leaf); */
    return leaf;
  }

  int_vertex++;
  this_vertex = int_vertex;
  k=0;
  for(i=1;i<len-1;i++) {
    if(s[i]=='(') k++;
    if(s[i]==')') k--;
    if(k==0) break;
  }
  j=i;
  

  s1 = (char*) space(sizeof(char)*(j+1));
  for(i=1;i<=j;i++) s1[i-1]=s[i];
  s1[j]='\0';
  
  s2 = (char *) space(sizeof(char)*(len-1-j));
  for(i=j+1;i<len-1;i++) s2[i-(j+1)] = s[i];
  s2[i-(j+1)]='\0';

  /* printf("> %s = %s + %s at <%d>\n",s,s1,s2,int_vertex); */

  e1 = fill_T(T,s1,0);
  T->E[e1].v2 = -this_vertex;
  e2 = fill_T(T,s2,0);
  T->E[e2].v2 = -this_vertex;
  
  free(s1);
  free(s2);

  int_edge++;

  T->I[this_vertex].e1 = e1;
  T->I[this_vertex].e2 = e2;
  T->I[this_vertex].e3 = int_edge;
  T->E[int_edge].v1 = -this_vertex;

  return int_edge;
  
  /* interior vertex */
  
}
char *InteriorTreeString(char *s)
{
  int i,j,n,N;
  char *si;
  n= strlen(s);
  N= number_of_leaves(s);
  si = (char*) space(sizeof(char)*(2*(N-2)+1));
  j=0;
  for(i=0;i<n;i++) {
    if( (s[i]=='(')&&(s[i+1]=='(')) { si[j] = '('; j++; }
    if( (s[i]==')')&&(s[i-1]==')')) { si[j] = ')'; j++; }
  }
  si[j]='\0';
  return si;
}
	
      
