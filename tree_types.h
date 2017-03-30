/* trees.h */

typedef struct {
  int e1;
  int e2;
  int e3;
} Node;

typedef struct {
  int v1;
  int v2;
} Edge;

typedef struct {
  int size;
  Node *I;
  Node *L;
  Edge *E;
} Tree;


/* End of file */
