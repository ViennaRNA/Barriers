/*
  Last changed Time-stamp: <2006-07-24 17:15:20 xtof>
  $Id: poho_generate.c,v 1.1 2006/07/25 14:34:55 xtof Exp $
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <sys/types.h>
#include <math.h>
#include "config.h"

#define  PI 3.1415927

static unsigned short xsubi[3];

static void nrerror(const char message[])       /* output message upon error */
{
  fprintf(stderr, "\n%s\n", message);
  exit(EXIT_FAILURE);
}

static void *space(unsigned size)
{
  void *pointer;

  if ( (pointer = (void *) calloc(1, (size_t) size)) == NULL) {
#ifdef EINVAL
    if (errno==EINVAL) {
      fprintf(stderr,"SPACE: requested size: %d\n", size);
      nrerror("SPACE allocation failure -> EINVAL");
    }
    if (errno==ENOMEM)
#endif
      nrerror("SPACE allocation failure -> no memory");
  }
  return  pointer;
}

#ifdef WITH_DMALLOC
#define space(S) calloc(1,(S))
#endif

static double urn(void)
/* uniform random number generator; urn() is in [0,1] */
/* uses a linear congruential library routine */
/* 48 bit arithmetic */
{
#ifdef HAVE_ERAND48
  extern double erand48(unsigned short[]);
  return erand48(xsubi);
#else
  return ((double) rand())/RAND_MAX;
#endif
}

void num2STRING(int dec, int l, char s[]);

int main(int argc, char *argv[])
{

  int i,j,k,n,nP,nF,l,N;
  char *s;
  char **P;
  int ***J;
  int seedval;
  int *EN;

  xsubi[0] = 4711;
  xsubi[1] = 42;
  xsubi[2] = 333;
  
  n = 10;   /* string length */
  nP = 3;   /* number of patterns to superimpose */
  nF = 2;   /* number of objective functions */
  
  if(argc>1)  n = atoi(argv[1]);        /* length */
  if(argc>2)  nF = atoi(argv[2]);       /* functions */
  if(argc>3)  nP = atoi(argv[3]);       /* patterns */
  if(argc>4)  {
    seedval = atoi(argv[4]);
    xsubi[0] = seedval%3417;
    xsubi[1] = seedval%4711;
    xsubi[2] = seedval%2347;
  }

  if(nP<=nF) {
    fprintf(stderr, "Warning, nP must be at least nF+1\n");
    nP=nF+1;
  }
  
  N = (int) pow(2,n);

  /* Setup */ 

  s = (char *) space((n+1)*sizeof(char));

  P = (char **) space((nP+nF-1)*sizeof(char*));
  for(i=0;i<nP+nF-1;i++) {
    P[i]= (char *) space(n*sizeof(char*));
    for(j=0;j<n;j++) {
      if( urn()>0.5 ) P[i][j] = '-';
      else P[i][j] = '+';
    }
    fprintf(stderr,"%s\n",P[i]);
  }

  J = (int ***) space(sizeof(int **)*nF);
  for(l=0;l<nF;l++) { 
    J[l] = (int **) space(sizeof(int *)*n);
    for(i=0;i<n;i++) J[l][i] = (int *) space(sizeof(int)*n);
    
    for(k=0;k<nP;k++) { 
      for(i=0;i<n;i++) {
	for(j=0;j<n;j++) {
	  if(P[k+l][i]==P[k+l][j])  J[l][i][j]++; 
	  else                      J[l][i][j]--;
	}
      }
    }
  }
  for(i=0;i<nP+nF-1;i++) free(P[i]);
  free(P);
  
  EN = (int *) space(sizeof(int)*nF);
  for(i=0;i<N;i++) {
    int energy,ss;
    num2STRING(i,n,s);
    printf("%s ",s);
    energy = 0;
    for(l=0;l<nF;l++) { 
      EN[l] =0;
      for(k=0; k<n;k++) {
	for(j=0;j<n;j++) {
	  ss = 1;	  
	  if(s[k]=='-') ss=-ss;
	  if(s[j]=='-') ss=-ss;
	  EN[l] -= J[l][k][j]*ss;
	}
      }
      energy += EN[l];
    }
    printf("%9.6f", (double)energy/(double)(n*nP*nF));
    for(l=0;l<nF;l++) printf(" %6d",EN[l]);
    printf("\n");  
  }
  free(EN);
  free(s);

  for(i=0;i<n;i++) printf("%c",'@');
  printf(" -99999999999999999  "); 
  printf(" P:%d Q2 H-%d \n", nF, nP);  

  for(l=0;l<nF;l++) {
     for(i=0;i<n;i++) free(J[l][i]);
     free(J[l]);
  }
  free(J);

  return 0;
} 

void num2STRING(int dec, int l, char s[])
{
   int i,d;
   char RNAlphaBet[] = "+-";

   d = dec ;
   for (i = l-1; i >= 0; i--) {
    
    s[i] = RNAlphaBet[d%2];

    d = (int) (d/2);
  }
  s[l]='\0';
}

/* *********************************************************************** */




