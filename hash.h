#define HASHSIZE 1048576   /* 2^20 */ 

void *hashtab[HASHSIZE+1];

extern int hash_comp(void *x, void *y);

