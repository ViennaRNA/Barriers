/***************************************************************************
          constraints.cpp  -  functions dealing with sequence constraints
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/

#include "secis_constraints.h"

/*********************************************************
Test whether a sequence consists only of valid symbols
concerning the IUPAC-code
*********************************************************/

int Check_iu(char *secis_nuc)
{
  int secis_size;
  secis_size = strlen(secis_nuc);
  
   for(int i=0; i<secis_size; i++)
   {
      secis_nuc[i] = toupper(secis_nuc[i]);
      if  (Valid_IUPAC(secis_nuc[i]) == 0)
         return 0;
   }
   return 1;
}



/*********************************************************
*   Test whether secis_nuc and secis_struct harmonize    *
*********************************************************/

int Check_secis_nuc_and_struct(char *secis_nuc, int *base_pairs)
{
  int bp_pos_j, bp_assign_i, bp_assign_j, secis_size;
  bool possible;

  secis_size = strlen(secis_nuc);
  
  for (int bp_pos_i=0; bp_pos_i<secis_size; bp_pos_i++) {
    bp_pos_j = base_pairs[bp_pos_i];

    if (bp_pos_j != -1) {  // if there is a base pair
      possible = false;
      for (int bp_assign = 0; bp_assign<6; bp_assign++) {
	BP2_2(bp_assign, bp_assign_i, bp_assign_j);
	
	if (Compare_IUPAC_Base(secis_nuc[bp_pos_i],bp_assign_i)
	    + Compare_IUPAC_Base(secis_nuc[bp_pos_j],bp_assign_j) == 2) {
	  possible = true;
	  break;
	}
      }
      
      if (possible == false) {
	cerr << "\n"    << secis_nuc[bp_pos_i]
	     << " and " << secis_nuc[bp_pos_j]
	     << " are not compatible!\n\n";
	exit(1);
      }
    }
  }
  
  return 1;
}



/******************************************************************
* Translation of the IUPAC-Seq (secis_nuc) to a twodimensional    *
* array, i.e. R becomes 1010, i.e. A and G are valid, C and U not *
******************************************************************/

int** getSeqConstraints(char* secis_nuc)
{
  int secis_size;
  int** seq_constraints;

  secis_size = strlen(secis_nuc);
  
  // allocate and initialize
  seq_constraints = (int**) malloc(sizeof(int*)*secis_size);
  for (int i=0; i<secis_size; i++) {
    seq_constraints[i] = (int*) malloc(sizeof(int)*4);
    for (int j=0; j<4; j++)
      seq_constraints[i][j] = Compare_IUPAC_Base(secis_nuc[i], j);
  }

  return (seq_constraints);
}
