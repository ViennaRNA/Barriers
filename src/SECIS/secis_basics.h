/***************************************************************************
                       basics.h  -  global variables and functions
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
 ***************************************************************************/

#ifndef _BASICS__
#define _BASICS__

#include <ctype.h>
#include <string.h>
#include <iostream>
#include <string>
#include <math.h>

using namespace std;

/******************************************************************************
*                                   Macros                                    *
******************************************************************************/
#define NUM_ACIDS 20
#define NUM_PAM 21


/******************************************************************************
*                               Constants                                     *
******************************************************************************/

extern const int aa_code[4][4][4];
extern const char int2char_aa[NUM_PAM];
extern const int BLOSUM[NUM_PAM][NUM_PAM];
extern const int PAM[NUM_PAM][NUM_PAM];

/******************************************************************************
*                         Functions                                           *
******************************************************************************/

int char2int_aa(char aa);
int* char2int_aa_seq(char* aa);
int sim_blosum(int* seq1, int* seq2, int len1, int len2);
int sim_blosum(char* seq1, char* seq2);
int aa_equal(int aa1, int aa2);
int aa_equal(char aa1, char aa2);
int count_aa_mutations(int* seq1, int* seq2, int len1, int len2);
int count_aa_mutations(char* seq1, char* seq2);

int* nucint_to_aaint(int* int_nuc_seq, int len);
int* nucchar_to_aaint(char* nuc_seq);

char int2char(int base);
char* int2char(int* num_seq, int size);
int char2int(char base);
int* char2int(char* seq);
void BP2_2(int bp, int &base_i, int &base_j);
int BP2int(int b1, int b2);
void BP2_2char(int bp, char &base_i, char &base_j);
int BPchar2int(char b1, char b2);

int Valid_IUPAC(char iu);
int Compare_IUPAC_Base(char iu, int base);
int Sum_INT_MIN2(int s1, int s2);
void print_int_seq(int* seq, int len);

#endif   // _BASICS_
