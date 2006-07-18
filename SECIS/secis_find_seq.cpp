/*
  Last changed Time-stamp: <2006-07-17 11:56:35 xtof>
  $Id: secis_find_seq.cpp,v 1.1 2006/07/18 14:08:49 xtof Exp $
*/
/***************************************************************************
                        find_seq.cpp  -  main
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
 ***************************************************************************/

#include <iostream>
#include <stdlib.h>

#include "secis_constraints.h"
#include "secis_generate.h"
#include "secis_struct.h"
extern "C" {
#include "utils.h"
}
  
extern long number;

using namespace std;

/******************************************************************************
*                                help for usage                               *
******************************************************************************/

void usage(char* name)
{
   cout << "\ncall: " << name
	<< " \"SECISstructure\" \"SECIS_IUPAC\" \"protein sequence\"\n";
   cout << "              [-M max. number of accepted amino acid mutations]\n";
   cout << "              [-S min. alignment score (amino acid sequence)]\n\n";
   exit(1);
}


/******************************************************************************
*                         Printing the input-values                           *
******************************************************************************/

void print_parameters(char* secis_nuc, char* secis_struct, char* protein_seq,
		      int max_mut, int min_algn_score)
{
   cout << endl;
   cout << "*******************************************" << endl;
   cout << "Parameters:" << endl;
   cout << "*******************************************" << endl;
   cout << "Structure  : " << secis_struct << endl;
   cout << "Constraints: " << secis_nuc << endl << endl;
   cout << "Protein sequence: " << protein_seq << endl << endl;
   cout << "Max. number of mutations : " << max_mut << endl;
   cout << "Min. alignment score (aa): " << min_algn_score << endl;
   cout << "*******************************************" << endl << endl;
}


/******************************************************************************
*                                 main-function                               *
******************************************************************************/

int main(int argc, char *argv[])
{
  char* secis_struct;    // structure in bracket notation and G-U for GU-pairs
  char* secis_nuc;       // IUPAC information about SECIS element
  char* protein_seq;     // amino acid sequence of the protein
  int* int_seq;
  int max_mut;           // max. number of accepted mutations in the
		         // amino acid sequence
  int min_algn_score;    // min. valid alignment score (amino acid)
  int protein_size;      // size of the given aa sequence
  int secis_size;

  
  max_mut = INT_MAX;
  min_algn_score = INT_MIN;

  /*
  if (argc < 3) {
    usage(argv[0]);
    exit(-1);
  }


  secis_struct = argv[1];
  secis_nuc    = argv[2];
  protein_seq  = argv[3];
  */
  secis_nuc = get_line(stdin);
  secis_struct = get_line(stdin);
  protein_seq = get_line(stdin);

  // process command line options
  for (int i = 1; i<argc; i++) {
    if (argv[i][0] == '-')
      switch (argv[i][1]) {
      case 'M': i++;
	if (argv[i] == NULL)
	  usage(argv[0]);
	max_mut = atoi(argv[i]);
	break;
      case 'S': i++;
	if (argv[i] == NULL)
	  usage(argv[0]);
	min_algn_score = atoi(argv[i]);
	break;
      default : usage(argv[0]);
      }
  }

  /****************************************************************
  *****************************************************************
  * checking for errors                                           *
  *****************************************************************
  ****************************************************************/

  // canonize input data
  secis_size = (int)strlen(secis_struct);
  for (int i=0; i<secis_size; i++)
    secis_struct[i] = toupper(secis_struct[i]);

  protein_size = (int)strlen(protein_seq);
  for (int i=0; i<protein_size; i++)
    protein_seq[i] = toupper(protein_seq[i]);

  // test whether secis_struct and secis_nuc are of the same length
  //**************************************************************
  if (secis_size != (int)strlen(secis_nuc)) {
    cerr << "\nThe structure and the nucleotide sequence of the "
	 << "SECIS element must have the same length!\n\n";
    exit(1);
  }
  
  // test whether structure is valid
  //*********************************
  if (! check_brackets(secis_struct)) {
    cerr << "\nNo valid structure!\n\n";
    exit(1);
  }
  
  // test whether secis_nuc contains valid symbols, if so, translate
  // IUPAC-Code to an array of valid bases
  //**************************************************************************
  if (! Check_iu(secis_nuc)) {
    cerr << "\nNo valid IUPAC Code!\n\n";
    exit(1);
  }

  // test whethe the protein sequence is at least as long as the
  // translated nucleotide sequence
  //**************************************************************************
  if (secis_size > 3*protein_size) {
    cerr << "\nThe given protein sequence is too short!\n\n";
    exit(1);
  }

  // test whether the max. number of accepted amino acid mutations ist valid
  //***********************************************************************
  if (max_mut < 0) {
    cerr << "\nThe number of accepted amino acid mutations has "
	 << "to be higher or equal to 0!\n";
    exit(1);
  }

   /****************************************************************
   *****************************************************************
   * The programm itself                                           *
   *****************************************************************
   ****************************************************************/

   // algorithm:
   //************
   print_parameters(secis_nuc, secis_struct, protein_seq,
		    max_mut, min_algn_score);

   int_seq = init_recursion(secis_nuc, secis_struct, protein_seq,
			    max_mut, min_algn_score);
   
   gen_seq(0,int_seq);

   printf("number of sequences: %ld\n", number);

   return (EXIT_SUCCESS);
}

