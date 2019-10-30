/***************************************************************************
                        basics.cpp  -  global variables and functions
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/
#include <climits>
#include "secis_basics.h"


/*****************************************************************
*    3D array, giving for each three nuc (int) an int            *
*    aa_code[codon[1]][codon[2]][codon[3]], zB: CGU=4            *
*****************************************************************/

const int aa_code[4][4][4] =
{
{
{11,2,11,2},
{16,16,16,16},
{1,15,1,15},
{9,9,12,9}
},
{
{5,8,5,8},
{14,14,14,14},
{1,1,1,1},
{10,10,10,10}
},
{
{6,3,6,3},
{0,0,0,0},
{7,7,7,7},
{19,19,19,19}
},
{
{20,18,20,18},
{15,15,15,15},
{20,4,17,4},
{10,13,10,13}
}
};
// 20 represents Stop-codons

/*****************************************************************
*    amino acids in one character representation                 *
*****************************************************************/

const char int2char_aa[NUM_PAM] =
{
   'A', 'R', 'N', 'D',
   'C', 'Q', 'E', 'G',
   'H', 'I', 'L', 'K',
   'M', 'F', 'P', 'S',
   'T', 'W', 'Y', 'V',
   '#'
};

/*****************************************************************
*                Similarity-Matrix (BLOSUM62)                    *
*   (plus a row and a column for the Stop-Sign, this is not      *
*   allowed)                                                     *
*****************************************************************/

const int BLOSUM[NUM_PAM][NUM_PAM] =
{
{4, -1, -2, -2, 0, -1, -1, 0, -2, -1, -1, -1, -1, -2, -1, 1, 0, -3, -2, 0, INT_MIN},
{-1, 5, 0, -2, -3, 1, 0, -2, 0, -3, -2, 2, -1, -3, -2, -1, -1, -3, -2, -3, INT_MIN},
{-2, 0, 6, 1, -3, 0, 0, 0, 1, -3, -3, 0, -2, -3, -2, 1, 0, -4, -2, -3, INT_MIN},
{-2, -2, 1, 6, -3, 0, 2, -1, -1, -3, -4, -1, -3, -3, -1, 0, -1, -4, -3, -3, INT_MIN},
{0, -3, -3, -3, 9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, INT_MIN},
{-1, 1, 0, 0, -3, 5, 2, -2, 0, -3, -2, 1, 0, -3, -1, 0, -1, -2, -1, -2, INT_MIN},
{-1, 0, 0, 2, -4, 2, 5, -2, 0, -3, -3, 1, -2, -3, -1, 0, -1, -3, -2, -2, INT_MIN},
{0, -2, 0, -1, -3, -2, -2, 6, -2, -4, -4, -2, -3, -3, -2, 0, -2, -2, -3, -3, INT_MIN},
{-2, 0, 1, -1, -3, 0, 0, -2, 8, -3, -3, -1, -2, -1, -2, -1, -2, -2, 2, -3, INT_MIN},
{-1, -3, -3, -3, -1, -3, -3, -4, -3, 4, 2, -3, 1, 0, -3, -2, -1, -3, -1, 3, INT_MIN},
{-1, -2, -3, -4, -1, -2, -3, -4, -3, 2, 4, -2, 2, 0, -3, -2, -1, -2, -1, 1, INT_MIN},
{-1, 2, 0, -1, -3, 1, 1, -2, -1, -3, -2, 5, -1, -3, -1, 0, -1, -3, -2, -2, INT_MIN},
{-1, -1, -2, -3, -1, 0, -2, -3, -2, 1, 2, -1, 5, 0, -2, -1, -1, -1, -1, 1, INT_MIN},
{-2, -3, -3, -3, -2, -3, -3, -3, -1, 0, 0, -3, 0, 6, -4, -2, -2, 1, 3, -1, INT_MIN},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4, 7, -1, -1, -4, -3, -2, INT_MIN},
{1, -1, 1, 0, -1, 0, 0, 0, -1, -2, -2, 0, -1, -2, -1, 4, 1, -3, -2, -2, INT_MIN},
{0, -1, 0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1, 1, 5, -2, -2, 0, INT_MIN},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1, 1, -4, -3, -2, 11, 2, -3, INT_MIN},
{-2, -2, -2, -3, -2, -1, -2, -3, 2, -1, -1, -2, -1, 3, -3, -2, -2, 2, 7, -1, INT_MIN},
{0, -3, -3, -3, -1, -2, -2, -3, -3, 3, 1, -2, 1, -1, -2, -2, 0, -3, -1, 4, INT_MIN},
{INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN}
};


/*****************************************************************
*                Similarity-Matrix (PAM250)                      *
*   (plus a row and a column for the Stop-Sign, this is not      *
*   allowed)                                                     *
*****************************************************************/

const int PAM[NUM_PAM][NUM_PAM] =
{
{2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -4, 1, 1, 1, -6, -3, 0, INT_MIN},
{-2, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, INT_MIN},
{0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -4, -1, 1, 0, -4, -2, -2, INT_MIN},
{0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, INT_MIN},
{-2, -4, -4, -5, 4, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, INT_MIN},
{0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, INT_MIN},
{0, -1, 1, 3, -5, 2, 4, 0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, INT_MIN},
{1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3, -5, -1, 1, 0, -7, -5, -1, INT_MIN},
{-1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0, -2, INT_MIN},
{-1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, INT_MIN},
{-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, INT_MIN},
{-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, INT_MIN},
{-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6, 0, -2, -2, -1, -4, -2, 2, INT_MIN},
{-4, -4, -4, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5, -3, -2, 0, 7, -1, INT_MIN},
{1, 0, -1, -1, -3, 0, -1, -1, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, INT_MIN},
{1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 3, 1, -2, -3, -1, INT_MIN},
{1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -2, 0, 1, 3, -5, -3, 0, INT_MIN},
{-6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, INT_MIN},
{-3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10, -2, INT_MIN},
{0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4, INT_MIN},
{INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN, INT_MIN},
};


/*****************************************************************
* returns the int of an amino acid given in the one-char-repres  *
*****************************************************************/

int char2int_aa(char aa)
{
   char tmp = toupper(aa);
   for (int i=0; i<NUM_PAM; i++)
      if (tmp == int2char_aa[i])
         return i;
   return INT_MAX;
}


/*****************************************************************
* returns the int-seq of an amino acid sequence given in the     *
* one-char-repres                                                *
*****************************************************************/

int* char2int_aa_seq(char* aa)
{
   int* int_seq = (int*) malloc(sizeof(int)*((int)strlen(aa)));
   for (int i=0; i<(int)strlen(aa); i++)
      int_seq[i] = char2int_aa(aa[i]);
   return int_seq;
}


/******************************************************************
* returns the similarity of 2 amino acid sequences (given as int) *
* (without gaps)                                                  *
******************************************************************/

int sim_blosum(int* seq1, int* seq2, int len1, int len2)
{
   int len = (len1 < len2 ? len1 : len2);
   int sim = 0;
   for (int i=0; i<len; i++)
      sim = Sum_INT_MIN2(sim,BLOSUM[seq1[i]][seq2[i]]);
   return sim;
}


int sim_blosum(char* seq1, char* seq2)
{
   int len = (strlen(seq1) < strlen(seq2) ? (int)strlen(seq1) : (int)strlen(seq2));
   int sim = 0;
   for (int i=0; i<len; i++)
      sim = Sum_INT_MIN2(sim,BLOSUM[char2int(seq1[i])][char2int(seq2[i])]);
   return sim;
}



/*****************************************************************
* tests whether 2 amino acids are equal (int)                    *
*****************************************************************/

int aa_equal(int aa1, int aa2)
{
   if (aa1 == aa2)
      return 1;
   else
      return 0;
}

int aa_equal(char aa1, char aa2)
{
   if (aa1 == aa2)
      return 1;
   else
      return 0;
}



/*****************************************************************
* counts the number of different aa-pos. in 2 aa sequences (int) *
* (if the sequences have diffent lengths, positions with no      *
*  counterpart are forgotten)                                    *
*****************************************************************/

int count_aa_mutations(int* seq1, int* seq2, int len1, int len2)
{
   int len = (len1 < len2 ? len1 : len2);
   int mut = 0;
   for (int i=0; i<len; i++)
      mut += (1-aa_equal(seq1[i],seq2[i]));
   return mut;
}


int count_aa_mutations(char* seq1, char* seq2)
{
   int len = (strlen(seq1) < strlen(seq2) ? (int)strlen(seq1) : (int)strlen(seq2));
   int mut = 0;
   for (int i=0; i<len; i++)
      mut += (1-aa_equal(seq1[i], seq2[i]));
   return mut;
}



/****************************************************************
*   translates an int-nucleotide sequence into an amino acid    *
*   sequence (int)                                              *
****************************************************************/

int* nucint_to_aaint(int* int_nuc_seq, int len)
{
   int aa_length = (int)(len/3);
   int* aa_seq = (int*) malloc(sizeof(int)*(aa_length+1));

   for (int i=0; i<aa_length; i++)
      aa_seq[i]
	= aa_code[int_nuc_seq[3*i]][int_nuc_seq[3*i+1]][int_nuc_seq[3*i+2]];

   return aa_seq;
}

/****************************************************************
*   translates an char-nucleotide sequence into an amino acid    *
*   sequence (int)                                              *
****************************************************************/

int* nucchar_to_aaint(char* nuc_seq)
{
   int nuc_length = (int)strlen(nuc_seq);
   int* int_nuc_seq = char2int(nuc_seq);

   return nucint_to_aaint(int_nuc_seq,nuc_length);
}


/**********************************************************
*  Translation of an integer_base to a character_base     *
**********************************************************/

char int2char(int base)
{
   if (base == 0) return 'A';
   else if (base == 1) return 'C';
   else if (base == 2) return 'G';
   else if (base == 3) return 'U';
   else
   {
      cerr << "There's something wrong in your int-base!";
      exit(1);
   }
}

/**********************************************************
Translation of a sequence of integer_bases to a sequence 
of characters
**********************************************************/

char* int2char(int* num_seq, int size)
{
  char* seq = (char*) calloc(size+1, sizeof(char));
  for (int i=0; i<size; i++)
    seq[i] = int2char(num_seq[i]);

  return (seq);
}


/**********************************************************
Translation of a character base to a integer base
**********************************************************/

int char2int(char base)
{
   if (toupper(base) == 'A') return 0;
   else if (toupper(base) == 'C') return 1;
   else if (toupper(base) == 'G') return 2;
   else if (toupper(base) == 'U') return 3;
   else
   {
      cerr << "There's something wrong in your char base!";
      exit(1);
   }
}


/**********************************************************
Translation of a sequence of character_bases to a sequence 
of integers
**********************************************************/

int* char2int(char* seq)
{
   int* num_seq = (int*) malloc(sizeof(int)*strlen(seq));
   for (int i=0; i<(int)strlen(seq); i++)
   {
      if (toupper(seq[i]) == 'A') num_seq[i]=0;
      else if (toupper(seq[i]) == 'C') num_seq[i]=1;
      else if (toupper(seq[i]) == 'G') num_seq[i]=2;
      else if (toupper(seq[i]) == 'U') num_seq[i]=3;
      else
      {
         cerr << "There's something wrong in your char-sequence!";
         exit(1);
      }
   }
   return num_seq;
}


/**********************************************************
Translation of a base pair given as an integer to
two integer base
**********************************************************/

void BP2_2(int bp, int &base_i, int &base_j)
{
   if (bp == 0) {base_i=0; base_j=3;}
   else if (bp == 1) {base_i=1; base_j=2;}
   else if (bp == 2) {base_i=2; base_j=1;}
   else if (bp == 3) {base_i=3; base_j=0;}
   else if (bp == 4) {base_i=2; base_j=3;}
   else if (bp == 5) {base_i=3; base_j=2;}
   else
   {
      cerr << "Wrong BasepairNr.!" << endl;
      exit(1);
   }
}


/**********************************************************
wandelt die 2 Basen eines BP in einen Integerwert um
**********************************************************/

int BP2int(int b1, int b2)
{
   int kind_of_BP = 0;
   if ( b1==0 && b2==3 ) kind_of_BP = 0;
   else if ( b1==1 && b2==2 ) kind_of_BP = 1;
   else if ( b1==2 && b2==1 ) kind_of_BP = 2;
   else if ( b1==3 && b2==0 ) kind_of_BP = 3;
   else if ( b1==2 && b2==3 ) kind_of_BP = 4;
   else if ( b1==3 && b2==2 ) kind_of_BP = 5;
   else
   {
      cerr << "Wrong Basepair!" << endl;
      exit(1);
   }
   return kind_of_BP;
}


/**********************************************************
Translation of a base pair given as an integer to
two character bases
**********************************************************/

void BP2_2char(int bp, char &base_i, char &base_j)
{
   if (bp == 0) {base_i='A'; base_j='U';}
   else if (bp == 1) {base_i='C'; base_j='G';}
   else if (bp == 2) {base_i='G'; base_j='C';}
   else if (bp == 3) {base_i='U'; base_j='A';}
   else if (bp == 4) {base_i='G'; base_j='U';}
   else if (bp == 5) {base_i='U'; base_j='G';}
   else
   {
      cerr << "Wrong BasepairNr.!" << endl;
      exit(1);
   }
}

/**********************************************************
wandelt die 2 Basen eines BP (given as char) in einen Integerwert um
**********************************************************/

int BPchar2int(char b1, char b2)
{
   int kind_of_BP = 0;
   if ( b1=='A' && b2=='U' ) kind_of_BP = 0;
   else if ( b1=='C' && b2=='G' ) kind_of_BP = 1;
   else if ( b1=='G' && b2=='C' ) kind_of_BP = 2;
   else if ( b1=='U' && b2=='A' ) kind_of_BP = 3;
   else if ( b1=='G' && b2=='U' ) kind_of_BP = 4;
   else if ( b1=='U' && b2=='G' ) kind_of_BP = 5;
   else
   {
      cerr << "Wrong Basepair!" << endl;
      exit(1);
   }
   return kind_of_BP;
}




/*************************************************************
Test whether a character is valid in the IUPAC Code
*************************************************************/
int Valid_IUPAC(char iu)
{
   if (   (iu == 'A') || (iu == 'C') || (iu == 'G') || (iu == 'U')
       || (iu == 'M') || (iu == 'R') || (iu == 'W') || (iu == 'S')
       || (iu == 'Y') || (iu == 'K') || (iu == 'V') || (iu == 'H')
       || (iu == 'D') || (iu == 'B') || (iu == 'N'))
      return 1;
   else
      return 0;
}

/*************************************************************
Test whether a given integer_base is valid for a given 
iupac_base
*************************************************************/
int Compare_IUPAC_Base(char iu, int base)
{
   if (base == 0)
   {
      if (   (iu == 'A') || (iu == 'M') || (iu == 'R') || (iu == 'W')
	  || (iu == 'V') || (iu == 'H') || (iu == 'D') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 1)
   {
      if (   (iu == 'C') || (iu == 'M') || (iu == 'S') || (iu == 'Y')
	  || (iu == 'V') || (iu == 'H') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 2)
   {
      if (   (iu == 'G') || (iu == 'R') || (iu == 'S') || (iu == 'K')
	  || (iu == 'V') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else if (base == 3)
   {
      if (   (iu == 'U') || (iu == 'W') || (iu == 'Y') || (iu == 'K')
	  || (iu == 'H') || (iu == 'D') || (iu == 'B') || (iu == 'N'))
         return 1;
      else
         return 0;
   }
   else
   {
      cerr << "Wrong Base!" << endl;
      exit(1);
   }
}


/*****************************************************************
* returns the sum of 2 integers, if one of them is MIN_INT,      *
* MIN_INT is returned                                            *
*****************************************************************/

int Sum_INT_MIN2(int s1, int s2)
{
   if ((s1 == INT_MIN) || (s2 == INT_MIN))
      return INT_MIN;
   else
      return (s1+s2);
}

/*****************************************************************
* prints an sequence of integers of length len                   *
*****************************************************************/

void print_int_seq(int* seq, int len)
{
   for(int i=0; i<len; i++)
      printf("%d",seq[i]);
   printf("\n");
}



