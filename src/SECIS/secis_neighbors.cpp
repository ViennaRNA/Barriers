#include <cstring>
#include <cstdlib>

#include "secis_constraints.h"
#include "secis_basics.h"
#include "secis_struct.h"
#include "secis_neighbors.h" // needed for correct name-mangling

extern "C" {
#ifdef _DEBUG_
#include "energy_const.h"
#include "fold_vars.h"
#include "fold.h"
#include "part_func.h"
#endif
#include "stapel.h" // neighbor stack from barriers
}

/* input data */
static char* secis_nuc;
static char* secis_struct;
static char* protein_seq;
static int max_mut;
static int min_algn_score;
static int* pair_table;
#ifdef _DEBUG_
static double  kT;
#endif

/* derived data */
static char* secis_brackets;
static int*  protein_seq_int;
static int** seq_constraints;
static int protein_size;


void initialize_SECIS(char* s_nuc, char* s_struct, char* p_seq,
		      int max_m, int min_as)
{
  int len;

  /*
   * memorize input data
   */
  
  /* secis nucleotide sequence */
  len = (int)strlen(s_nuc);
  secis_nuc = (char*)calloc(len+1, sizeof(char));
  strcpy(secis_nuc, s_nuc);
  
  /* secis secondary structure */
  len = (int)strlen(s_struct);
  secis_struct = (char *)calloc(len+1, sizeof(char*));
  strcpy(secis_struct, s_struct);

  /* wildtype protein sequence */
  len = (int)strlen(p_seq);
  protein_seq = (char *)calloc(len+1, sizeof(char*));
  strcpy(protein_seq, p_seq);
  
  /* make pair table */
  pair_table = make_BasePair_Table(secis_struct);
  
  /* maximal mutation distance from wildetype protein; minimum
     alignment score of the mutant protein with the wildtype
     protein */
  max_mut = max_m;
  min_algn_score = min_as;

#ifdef _DEBUG_
  kT = (temperature+K0)*GASCONST/1000.0;
#endif
  /*
   * calculate derived data
  */
  len = (int)strlen(secis_struct);
  secis_brackets = (char *)calloc(len+1, sizeof(char*));
  strcpy(secis_brackets, secis_struct);
  for (int i=0; i<len; i++)
  {
     if (toupper(secis_brackets[i]) == 'G')
        secis_brackets[i] = '(';
     else if (toupper(secis_brackets[i]) == 'U')
        secis_brackets[i] = ')';
  }

  seq_constraints = getSeqConstraints(secis_nuc);
  protein_seq_int = char2int_aa_seq(protein_seq);
  protein_size = (int)strlen(protein_seq);
  ini_stapel(len); // inizialize barriers stack
}

/**
 * the scooring function:
 * translated protein sequence:
 * (1) mutation distance from wildtype <= MAX_MUT
 * (2) similarity score to wildtype >= MIN_ALN_SCORE
 * RNA sequence:
 * no constraints yet
 */
int make_scores(char* rna_seq) //, int secis_size, int secis_aa_size)
{
  int secis_size, secis_aa_size;
  int* aa_seq_int;
#ifdef _DEBUG_
  float energy, ensemble_energy, prob;
  char *pf_struct;
#endif
  int aa_mut, algn_score, ok_toggle;

  secis_size = strlen(rna_seq);
  secis_aa_size = (int)secis_size/3;
  ok_toggle = 0;
   
  aa_seq_int = nucchar_to_aaint(rna_seq); //translate to aa-seq
  aa_mut = count_aa_mutations(aa_seq_int, protein_seq_int,
			      secis_aa_size, protein_size); //count mutations

  if (aa_mut <= max_mut) {
    algn_score = sim_blosum(aa_seq_int, protein_seq_int,
			    secis_aa_size, protein_size);  //alignment score
    if (algn_score >= min_algn_score) {
#ifdef _DEBUG_
      energy = energy_of_struct(rna_seq,secis_brackets);

      // partition function fold
      init_pf_fold(secis_size);
      pf_struct = (char*) malloc(sizeof(char)*(secis_size+1));
      ensemble_energy = pf_fold(rna_seq,pf_struct);
      free_pf_arrays();
      prob = exp(-(energy-ensemble_energy)/(kT));
      printf("%s %6.2f %.8f %3d %3d\n",
	     rna_seq, energy, prob, algn_score, aa_mut);
      free(pf_struct);
#endif
      ok_toggle = 1; /* neighbor meets our requirements */
    }
  }
  free(aa_seq_int);

  return (ok_toggle);
}

void SECIS_move_it(char* rna_seq)
{
   char *rna_seq_new;
   int secis_size, secis_aa_size, current_bp_assign, num_bp;

   secis_size = (int)strlen(rna_seq);
   secis_aa_size = (int)(secis_size/3);

   rna_seq_new = (char *)calloc(secis_size+1, sizeof(char*));
   strcpy(rna_seq_new, rna_seq);
   
   for (int pos=0; pos<secis_size; pos++) {
     if ((secis_struct[pos] == '(') || (toupper(secis_struct[pos]) == 'G')) {
       current_bp_assign = BPchar2int(rna_seq[pos], rna_seq[pair_table[pos]]);
       if (secis_struct[pos] == '(')
	 num_bp = 4;
       else
	 num_bp = 6;
         
       for (int bp_assign=0; bp_assign<num_bp; bp_assign++) {
	 if (bp_assign == current_bp_assign)
	   continue;
	 else {
	   /**
	    * rna_seq_new mutieren (kann vorher schon an dieser Stelle
	    * mutiert worden sein, aber es wird ja sowieso mit
	    * OriginalBP aus rna_seq verglichen (current_bp_assign)
	    */
	   BP2_2char(bp_assign, rna_seq_new[pos],
		     rna_seq_new[pair_table[pos]]);

	   // scores testen
	   //if (make_scores(rna_seq_new, secis_size, secis_aa_size))
	   if (make_scores(rna_seq_new))
	     push(rna_seq_new); // push sequence to barriers neighbor stack
	 }
       }
       /**
	* neue Seq. wieder auf Ausgangszustand setzen (d.h. identisch
	* mit rna_seq)
	*/
       BP2_2char(current_bp_assign, rna_seq_new[pos],
		 rna_seq_new[pair_table[pos]]);
     }
     else if (secis_struct[pos] == '.') {
       for (int base=0; base<4; base++) {
	 if (int2char(base) == rna_seq[pos])
	   continue;
	 else {
	   /**
	    * rna_seq_new mutieren (kann vorher schon an dieser Stelle
	    * mutiert worden sein, aber es wird ja sowieso mit
	    * OriginalBP aus rna_seq verglichen (current_bp_assign)
	    */
	   rna_seq_new[pos] = int2char(base);

	   // scores testen
	   // if (make_scores(rna_seq_new, secis_size, secis_aa_size))
	   if (make_scores(rna_seq_new))
	     push(rna_seq_new); // push sequence to barriers neighbor stack
	 }
       }
       
       /**
	* neue Seq. wieder auf Ausgangszustand setzen (d.h. identisch
	* mit rna_seq)
	*/
       rna_seq_new[pos] = rna_seq[pos];
     }
   }
}

#if 0
/**
 * C interface
 */
void
initialize_secis(char* s_nuc, char* s_struct, char* p_seq,
		 int max_m, int min_as)
{
  initialize_SECIS(s_nuc, s_struct, p_seq, max_m, min_as);
}

void
secis_move_it(char* rna_seq)
{
  SECIS_move_it(rna_seq);
}

#endif
/* End of file */
