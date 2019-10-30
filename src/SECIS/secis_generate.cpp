/***************************************************************************
                generate.cpp  -  generating valid sequences
                             -------------------
    begin                : 27-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/

#include "secis_constraints.h"
#include "secis_generate.h"
#include "secis_struct.h"

extern "C" {
#include "energy_const.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
}

/*********************************************************
* generates all valid sequences, calculates their mfe,   *
* their frequencies in the ensembles, their protein-     *
* alignment-scores, and the number of mutated amino acids*
* after translating the RNA-sequences                    *
*    -> prints all of that                               *
* and count the number of valid sequences                *
*********************************************************/


/******************************************************************************
*                               global variables                              *
******************************************************************************/
long number;

static char* secis_brackets;
static char* secis_struct;
static char* secis_nuc;
static char* protein_seq;
static int* int_seq;
static int* base_pairs;
static int* protein_seq_int;
static int** seq_constraints;
static int secis_size;
static int secis_aa_size;
static int protein_size;
static int max_mut;
static int min_algn_score;
static double  kT;


static void allocate_arrays(void) {
  secis_brackets = (char*) calloc(secis_size+1, sizeof(char));
  int_seq        = (int*)  calloc(secis_size,   sizeof(int));
  base_pairs = make_BasePair_Table(secis_struct);
  
  if (Check_secis_nuc_and_struct(secis_nuc, base_pairs))
    seq_constraints = getSeqConstraints(secis_nuc);

  protein_seq_int = char2int_aa_seq(protein_seq);  
}

int* init_recursion(char* s_nuc, char* s_struct, char* p_seq,
		    int max_m, int min_as)
{
  number = 0;
  secis_nuc = s_nuc;
  secis_struct = s_struct;
  protein_seq = p_seq;
  secis_size    = strlen(secis_nuc);
  secis_aa_size = (int)(secis_size/3);
  protein_size  = (int)strlen(protein_seq);
  kT = (temperature+K0)*GASCONST/1000.0;
  max_mut = max_m;
  min_algn_score = min_as;

  allocate_arrays();
  
  for (int i=0; i<secis_size; i++) {
    if (secis_struct[i] == 'G')
      secis_brackets[i] = '(';
    else if (secis_struct[i] == 'U')
      secis_brackets[i] = ')';
    else
      secis_brackets[i] = secis_struct[i];
  }

  
  // initialize start sequence with all A's
  for (int i=0; i<secis_size; i++) int_seq[i] = 0;

  return (int_seq);
}

void gen_seq(int pos, int* seq)
{
  int* aa_seq;
  float energy, ensemble_energy, prob;
  char *nuc_seq, *pf_struct;
  int aa_mut, algn_score, secis_size;

  secis_size = strlen(secis_struct);

  if ((secis_struct[pos] == '(') || (secis_struct[pos] == 'G')) {
    for (int nuc=0; nuc<4; nuc++) {
      seq[pos] = nuc;
      seq[base_pairs[pos]] = 3-nuc; //complementary base (no g-u so far)

      if (   (seq_constraints[pos][nuc] == 1)
	  && (seq_constraints[base_pairs[pos]][3-nuc] == 1)) {

	if (pos < secis_size-1) gen_seq(pos+1,seq);
	else {
	  aa_seq = nucint_to_aaint(seq, secis_size); //translate to aa-seq
	  // count mutations
	  aa_mut = count_aa_mutations(aa_seq, protein_seq_int,
				      secis_aa_size, protein_size);

	  if (aa_mut <= max_mut) {
	    // alignment score
	    algn_score = sim_blosum(aa_seq, protein_seq_int,
				    secis_aa_size, protein_size);

	    if (algn_score >= min_algn_score) {
	      nuc_seq = int2char(seq,secis_size); // int-RNA in char RNA
	      energy = energy_of_struct(nuc_seq,secis_brackets);

	      // partition function folding
	      init_pf_fold(secis_size);
	      pf_struct = (char*) calloc(secis_size+1, sizeof(char));
	      ensemble_energy = pf_fold(nuc_seq,pf_struct);
	      free_pf_arrays();
	      prob = exp(-(energy-ensemble_energy)/(kT));
	      // nuc_seq mfe prob aln_scoor aa_mut
	      printf("%s %6.2f %g %3d %3d\n",
		     nuc_seq, energy, prob, algn_score, aa_mut);

	      number++;
	      free(nuc_seq);
	      free(pf_struct);
	    }
	  }
	  free(aa_seq);
	}
      }
    }// end for (nuc)

    if (secis_struct[pos] == 'G') {
      for (int nuc=2; nuc<4; nuc++) {
	seq[pos] = nuc;
	seq[base_pairs[pos]] = 5-nuc; //complementary base in g-u-pair

	if (   (seq_constraints[pos][nuc] == 1)
	    && (seq_constraints[base_pairs[pos]][5-nuc] == 1)) {
	  
	  if (pos < secis_size-1) gen_seq(pos+1,seq);
	  else {
	    aa_seq = nucint_to_aaint(seq, secis_size); //translate to aa-seq
	    // count mutations
	    aa_mut = count_aa_mutations(aa_seq, protein_seq_int,
					secis_aa_size, protein_size);

	    if (aa_mut <= max_mut) {
	      // alignment score
	      algn_score = sim_blosum(aa_seq, protein_seq_int,
				      secis_aa_size, protein_size);

	      if (algn_score >= min_algn_score) {
		nuc_seq = int2char(seq,secis_size); // int-RNA in char RNA
		energy = energy_of_struct(nuc_seq,secis_brackets);

		// partition function folding
		init_pf_fold(secis_size);
		pf_struct = (char*) calloc(secis_size+1, sizeof(char));
		ensemble_energy = pf_fold(nuc_seq,pf_struct);
		free_pf_arrays();
		prob = exp(-(energy-ensemble_energy)/(kT));
		// nuc_seq mfe prob aln_scoor aa_mut
		printf("%s %6.2f %g %3d %3d\n",
		       nuc_seq, energy, prob, algn_score, aa_mut);

		number++;
		free(nuc_seq);
		free(pf_struct);
	      }
	    }
	    free(aa_seq);
	  }
	}
      }//end for (nuc)
    }//end if (== 'G')

  }// end if
  else if (secis_struct[pos] == '.') {
    for (int nuc=0; nuc<4; nuc++) {
      seq[pos] = nuc;

      if (seq_constraints[pos][nuc] == 1) {
	
	if (pos < secis_size-1) gen_seq(pos+1,seq);
	else {
	  aa_seq = nucint_to_aaint(seq, secis_size); //translate to aa-seq
	  // count mutations
	  aa_mut = count_aa_mutations(aa_seq, protein_seq_int,
				      secis_aa_size, protein_size);

	  if (aa_mut <= max_mut) {
	    // alignment score
	    algn_score = sim_blosum(aa_seq, protein_seq_int,
				    secis_aa_size, protein_size);

	    if (algn_score >= min_algn_score) {
	      nuc_seq = int2char(seq,secis_size); //int-RNA in char RNA
	      energy = energy_of_struct(nuc_seq,secis_brackets);

	      // partition function folding
	      init_pf_fold(secis_size);
	      pf_struct = (char*) calloc(secis_size+1, sizeof(char));
	      ensemble_energy = pf_fold(nuc_seq,pf_struct);
	      free_pf_arrays();
	      prob = exp(-(energy-ensemble_energy)/(kT));
	      // nuc_seq mfe prob aln_scoor aa_mut
	      printf("%s %6.2f %g %3d %3d\n",
		     nuc_seq, energy, prob, algn_score, aa_mut);

	      number++;
	      free(nuc_seq);
	      free(pf_struct);
	    }
	  }
	  free(aa_seq);
	}
      }
    }// end for (nuc)
  }
  else {
    
    if (pos < secis_size-1) gen_seq(pos+1,seq);
    else {
      aa_seq = nucint_to_aaint(seq, secis_size); //translate to aa-seq
      // count mutations
      aa_mut = count_aa_mutations(aa_seq, protein_seq_int,
				  secis_aa_size, protein_size);

      if (aa_mut <= max_mut) {
	// alignment score
	algn_score = sim_blosum(aa_seq, protein_seq_int,
				secis_aa_size, protein_size);

	if (algn_score >= min_algn_score) {
	  nuc_seq = int2char(seq,secis_size); //int-RNA in char RNA
	  energy = energy_of_struct(nuc_seq,secis_brackets);

	  // partition function folding
	  init_pf_fold(secis_size);
	  pf_struct = (char*) calloc(secis_size+1, sizeof(char));
	  ensemble_energy = pf_fold(nuc_seq,pf_struct);
	  free_pf_arrays();
	  prob = exp(-(energy-ensemble_energy)/(kT));
	  // nuc_seq mfe prob aln_scoor aa_mut
	  printf("%s %6.2f %g %3d %3d\n",
		 nuc_seq, energy, prob, algn_score, aa_mut);

	  number++;
	  free(nuc_seq);
	  free(pf_struct);
	}
      }
      free(aa_seq);
    }
  }
}

void check_consistancy_of_input_data(){
  secis_size = (int)strlen(secis_struct);
  secis_aa_size = (int)(secis_size/3);
  for (int i=0; i<secis_size; i++)
    secis_struct[i] = toupper(secis_struct[i]);

  secis_brackets = (char*) calloc(secis_size+1, sizeof(char));
  for (int i=0; i<secis_size; i++) {
    if (secis_struct[i] == 'G')
      secis_brackets[i] = '(';
    else if (secis_struct[i] == 'U')
      secis_brackets[i] = ')';
    else
      secis_brackets[i] = secis_struct[i];
  }

   protein_size = (int)strlen(protein_seq);
   for (int i=0; i<protein_size; i++)
      protein_seq[i] = toupper(protein_seq[i]);
   protein_seq_int = char2int_aa_seq(protein_seq);

   // test whether structure is valid
   //*********************************
   int correct = check_brackets(secis_struct);
   if (correct)
     base_pairs = make_BasePair_Table(secis_struct);
   else
   {
      cerr << "\nNo valid structure!\n\n";
      exit(1);
   }

   // test whether secis_struct and secis_nuc are of the same length
   //**************************************************************
   if (secis_size != (int)strlen(secis_nuc))
   {
      cerr << "\nThe structure and the nucleotide sequence of the "
	   << "SECIS element must have the same length!\n\n";
      exit(1);
   }

   // test whether secis_nuc contains valid symbols, if so, translate
   // IUPAC-Code to an array of valid bases
   //**************************************************************************
   if (Check_iu(secis_nuc))
   {
     if (Check_secis_nuc_and_struct(secis_nuc, base_pairs))
         seq_constraints = getSeqConstraints(secis_nuc);
   }
   else
   {
      cerr << "\nNo valid IUPAC Code!\n\n";
      exit(1);
   }

   // test whethe the protein sequence is at least as long as the
   // translated nucleotide sequence
   //**************************************************************************
   if (secis_size > 3*protein_size)
   {
      cerr << "\nThe given protein sequence is too short!\n\n";
      exit(1);
   }

   // test whether the max. number of accepted amino acid mutations ist valid
   //***********************************************************************
   if (max_mut < 0)
   {
      cerr << "\nThe number of accepted amino acid mutations has "
	   << "to be higher or equal to 0!\n";
      exit(1);
   }

}
