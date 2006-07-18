/*
  Last changed Time-stamp: <2006-07-17 11:35:26 xtof>
  $Id: secis_struct.cpp,v 1.1 2006/07/18 14:08:49 xtof Exp $
*/
/***************************************************************************
                    struct.cpp  -  functions dealing with the structure
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/
#include "secis_struct.h"

/*--------------------------------------------------------------------------*/
/************************************************
Test whether secis_struct is a valid structure
(the structure "..........." is not valid)

copied: Vienna RNA Package + extended:
************************************************/

int check_brackets(char *line)
{
  int i,o,bonds,gu;

  i=o=bonds=gu=0;
  while( line[i] ){
    switch(line[i]) {
    case '(' :
      o++;
      i++;
      bonds++;
      break;
    case 'G':
      gu++;
      i++;
      bonds++;
      break;
    case '.' :
      i++;
      break;
    case ')' : 
      i++;
      o--;
      if(o<0) return 0;
      break;
    case 'U' :
      i++;
      gu--;
      if (gu<0) return 0;
      break;
    default:
      return 0;
    }
  }
  if (o>0) return 0;
  if (gu>0) return 0;
  if (bonds == 0) return 0;
  return 1;
}

/*--------------------------------------------------------------------------*/

/*****************************************************
Finds for each base the pairing base, free bases 
are set to -1

copied: Vienna RNA Package + extended:
*****************************************************/

int* make_BasePair_Table(char *structure)
{
  int j,hx;
  int *stack, *table;

  hx=0;
  stack = (int *) malloc(sizeof(int)*(strlen(structure)+1));
  table = (int *) malloc(sizeof(int)*(strlen(structure)+1));

  for (unsigned int i=0; i<strlen(structure); i++) {
    switch (structure[i]) {
    case '.':
      table[i]= -1;
      break;
    case 'G':
    case '(':
      stack[hx++]=i;
      break;
    case 'U':
    case ')':
      j = stack[--hx];

      if (hx<0) {
	cerr << structure << endl;
	cerr << "unbalanced brackets in make_BasePair_Table";
      }

      if (   ((structure[i] == 'U') && (structure[j] != 'G'))
	  || ((structure[i] == ')') && (structure[j] != '('))) {
	cerr << structure << endl;
	cerr << "unbalanced G-U brackets in make_BasePair_Table";
	exit(1);
      }
      
      table[i]=j;
      table[j]=i;

      break;
    }
  }
  
  if (hx!=0) {
    cerr << structure << endl;
    cerr << "unbalanced brackets in make_BasePair_Table";
  }
  
  free(stack);

  return (table);
}




