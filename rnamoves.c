/*
 * Alternative move set definitions for RNAs
 */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "rnamoves.h"
#include "stapel.h"


extern int MYTURN;       /**< minimum loop length +1, excluding closing pair */
const  int UNPRD = -1;   /**< value assigned to unpaired bases in pairtable */
const  int VERB  =  1;   /**< be verbose? */


/** flag to track initialization state of RNA2_data_t structure object */
int isInitialized = 0;


/**
 * @brief data structure for global variables of move generator
 */
typedef struct RNA2_data {
    char *seq;      /**< sequence RNA sequence */
    int shift;      /**< shift flag, turn on generation of shift moves */
    int noLP;       /**< noLP flag, turn on move generation for noLP structures */

    unsigned int len; /**< length of sequence */

    char *form; /**< dot bracket structure for neighbors */

    /** pairing array for efficient traversal and checks
     *
     * Invariant:
     *  pt[i] =
     *    UNPRD (i.e. negative constant) if form[i]=='.';
     *  OR
     *    j if form[i]!='.' and i is paired to j
     *
     * The pairtable is zero-based, i.e. pt[0] is the pairing partner of first base
     */
    int *pt;        /**< pairtable, i.e. pt[i]==j iff base i pairs with base j

    /** local stack for structure parser */
    int *stack;
} RNA2_data_t;


/** global variables of move generator */
RNA2_data_t *r2d;


/**
 * @brief Check whether (a,b) is a valid basepair. Valid basepairs are the
 *  canonical pairs AU and GC and also the wobble pair GU (and their respective
 *  symmetrical counterparts UA, CG, UG).
 * @param a First base
 * @param b Second base
 * @return 1 if basepair is valid, 0 if not
 */
int valid_bp_char( char const a, char const b)
{
    if(
            a == 'A' && b == 'U' ||
            a == 'U' && b == 'A' ||
            a == 'G' && b == 'C' ||
            a == 'C' && b == 'G' ||
            a == 'G' && b == 'U' ||
            a == 'U' && b == 'G'
            )
        return 1;
    else
        return 0;
}


/**
 * @brief Check the sequence for validity, convert it to upper case and
 *  substitute Ts for Us. Exits if sequence contains invalid characters.
 */
void normalize_seq()
{
    int i;
    for( i=0; i < r2d->len; i++)
        switch( r2d->seq[i])
        {
        case 't':
        case 'T':
            r2d->seq[i] = 'U';
            break;
        case 'a':
        case 'u':
        case 'g':
        case 'c':
            r2d->seq[i] += 'A'-'a';    /* convert to upper case */
            break;
        case 'A':
        case 'U':
        case 'G':
        case 'C':
            break;                                  /* valid, do nothing */
        default:
            fprintf( stderr, "ERROR: normalize_seq(): invalid base '%c' at "
                             "position %d of sequence: \n%s\n",
                     r2d->seq[i], i, r2d->seq);
            exit(1);                                /* invalid, return error */
        }
    return;
}


/**
 * @brief Check whether (i,j) is a valid basepair in the stored sequence,
 *  c.f. valid_bp_char(). It is assumed that i and j are less than the sequence
 *  length.
 * @param i First index of stored sequence. i<r2d->len is assumed.
 * @param j Second index of stored sequence. j<r2d->len is assumed.
 * @return 1 if basepair is valid, 0 if not
 */
inline int valid_bp( int const i, int const j)
{
    return valid_bp_char( r2d->seq[i], r2d->seq[j]);
}


void RNA2_init(char const * const sequence, int const shift, int const noLP)
{
    int i;
    if ( isInitialized)     /* do not leak anything if already intialized */
        RNA2_free();
    isInitialized = 1;
    r2d = (RNA2_data_t *) malloc( sizeof(RNA2_data_t));

    r2d->len = strlen(sequence);
    r2d->seq = (char*) malloc( (r2d->len+1)*sizeof(char));
    strcpy( r2d->seq, sequence);
    normalize_seq();   /* make uppercase, convert T to U */

    r2d->shift=shift;
    r2d->noLP=noLP;

    r2d->form = (char *)malloc((r2d->len+1) * sizeof(char));

    r2d->pt = (unsigned int *)malloc((r2d->len+1)*sizeof(unsigned int));
    r2d->stack = (unsigned int*)malloc(r2d->len*sizeof(unsigned int));

    for( i=0; i<r2d->len; i++) /* initialize with open structure */
    {
        r2d->form[i] = '.';
        r2d->pt[i] = UNPRD;
    }
    r2d->form[r2d->len] = 0; /* set terminating 0 */

    return;
}


/**
  @brief write parenthesis into form, adjust pairtable. Well-formedness of
   basepair (i,j) is assumed (non-crossing, valid basepair).
  @param i left end
  @param i right end
*/
inline void close_bp(const unsigned int i, const unsigned int j) {
    r2d->form[i]='(';
    r2d->form[j]=')';
    r2d->pt[i] = j;
    r2d->pt[j] = i;
    return;
}


/**
  @brief make base pair ends unpaired in form, adjust pairtable
  @param i left end
  @param i right end
*/
inline void open_bp(const unsigned int i, const unsigned int j) {
    r2d->form[i] = r2d->form[j] = '.';
    r2d->pt[i] = r2d->pt[j] = UNPRD;
    return;
}


inline int seq_len()
{
    return r2d->len;
}


inline char* dot_bracket_str()
{
    return r2d->form;
}


inline int is_paired(int i)
{
    return r2d->form[i] != '.';
}


inline int is_opening_bp(int i)
{
    return r2d->form[i] == '(';
}


inline int is_closing_bp(int i)
{
    return r2d->form[i] == ')';
}


inline int pairs_with(int i)
{
    return r2d->pt[i];
}


inline int is_bp(int i, int j)
{
    return r2d->pt[i] == j;
}


inline int has_min_loop_length (int i, int j)
{
    return j-i >= MYTURN;
}


/**
 * @brief Checks whether basepair (i,j) is outside-lonely, i.e. whether
 *  (i-1,j+1) is not a basepair in the structure.
 * @param i left index
 * @param j right index
 * @return 1 if (i,j) is outside-lonely, 0 if not
 */
inline int is_out_lonely( int i, int j)
{
    return i==0 || !is_bp(i-1, j+1);
}


/**
 * @brief Checks whether basepair (i,j) is inside-lonely, i.e. whether
 *  (i+1,j-1) is not a basepair in the structure.
 * @param i left index
 * @param j right index
 * @return 1 if (i,j) is inside-lonely, 0 if not
 */
inline int is_ins_lonely( int i, int j)
{
    return !is_bp(i+1, j-1);
}


/**
 * @brief Checks whether basepair (i,j) is lonely, i.e. whether it is inside-
 *  or outside-lonely.
 * @param i left index
 * @param j right index
 * @return 1 if (i,j) is lonely, 0 if not
 */
inline int isLonely( int i, int j)
{
    return is_out_lonely( i, j) && is_ins_lonely( i, j);
}


/**
 * @brief Checks whether (i+1,j-1) is a basepair that would grow lonely (i.e.
 *  becomes a lonely pair) if the basepair (i,j) was removed. It is not required
 *  that (i,j) is currently contained in the structure.
 * @param i left index
 * @param j right index
 * @return 1 if (i+1,j-1) is a basepair in the current structure and grows lonely
 *  when removing (i,j), 0 if not
 */
inline int insGrowLonely( int i, int j)
{
    return is_bp(i+1, j-1) && is_ins_lonely(i+1, j-1);
}


/**
 * @brief Checks whether (i-1,j+1) is a basepair that would grow lonely (i.e.
 *  becomes a lonely pair) if the basepair (i,j) was removed. It is not required
 *  that (i,j) is currently contained in the structure.
 * @param i left index
 * @param j right index
 * @return 1 if (i-1,j+1) is a basepair in the current structure and grows lonely
 *  when removing (i,j), 0 if not
 */
inline int outGrowLonely( int i, int j)
{
    return i!=0 && is_bp(i-1, j+1) && is_out_lonely(i-1, j+1);
}


/**
 * @brief Given an unpaired index i, find the first index j with i<j such that (i,j)
 *  is a valid basepair in the current sequence and structure. To get all such
 *  positions j, repeatedly call this function passing the index j returned by
 *  the last call.
 * @param i Left index, assumed to be unpaired
 * @param j Last index such that (i,j) is a valid basepair. To get the first
 *  such j, pass j=i. i<=j is assumed.
 * @return next index k with i<=j<k such that (i,k) is a valid basepair, or
 *  UNPRD if no such k exists
 */
int RNA2_move_nxt_val_bp_right( int i, int j)
{
    do /* jump over inner base pairs of current loop */
    {
        j++;
        while(j<seq_len() && is_opening_bp(j))
            j = pairs_with(j)+1;
        if(j>=seq_len() || is_closing_bp(j))  /* here, the current loop ends */
            return UNPRD;
    } while( !has_min_loop_length(i, j) || !valid_bp(i,j));

    return j;
}


/**
 * @brief Generate neighbors by inserting every possible basepair. Only canonical
 *  structures (i.e. ones not containing lonely pairs) are generated. If
 *  necessary, lonely 2-stacks are inserted. The resulting structures are pushed
 *  onto the global structure stack.
 */
void RNA2_move_noLP_bpins() {
    int i,j;
    for( i=0; i<seq_len(); i++)
        if( !is_paired(i))
        {
            j=i;
            while( (j = RNA2_move_nxt_val_bp_right(i, j)) != UNPRD)
            {
                close_bp( i, j);
                if( isLonely( i, j))
                {
                    if( has_min_loop_length(i+1, j-1)
                        && !is_paired(i+1)
                        && !is_paired(j-1)
                        && valid_bp(i+1, j-1)
                        && is_ins_lonely(i+1, j-1)
                    )
                    {   /* insert lonely stack */
                        close_bp( i+1, j-1);
                        if(VERB)
                            fprintf( stderr, "pushing lsi %s\n", dot_bracket_str());
                        push(dot_bracket_str());
                        open_bp( i+1, j-1);
                    }
                }
                else    /* bp not lonely, insert*/
                {
                    if(VERB)
                        fprintf( stderr, "pushing lpi %s\n", dot_bracket_str());
                    push(dot_bracket_str());
                }
                open_bp( i, j);
            }
        }
    return;
}


/**
 * @brief Generate neighbors by deleting every possible basepair. Only canonical
 *  structures (i.e. ones not containing lonely pairs) are generated. If
 *  necessary, lonely 2-stacks are removed. The resulting structures are pushed
 *  onto the global structure stack.
 */
void RNA2_move_noLP_bpdel() {
    int i,j, ilg;
    for( i=0; i<seq_len(); i++)
        if(is_opening_bp(i))
        {
            j = pairs_with(i);
            open_bp( i, j);
            ilg = insGrowLonely( i, j);     /* is structure inside-lonely-growing? */
            if( ilg || outGrowLonely( i,j))
            {
                if( ilg && is_out_lonely(i, j) ) /* only delete lonely stack on the inside */
                {
                    open_bp( i+1, j-1);
                    if(VERB)
                        fprintf( stderr, "pushing lsd %s\n", dot_bracket_str());
                    push(dot_bracket_str());
                    close_bp( i+1, j-1);
                }
            }
            else    /* no lonely pairs arise */
            {
                if(VERB)
                    fprintf( stderr, "pushing lpd %s\n", dot_bracket_str());
                push(dot_bracket_str());
            }
            close_bp( i, j);
        }
    return;
}


/**
 * @brief Generate neighbors by shifting every possible basepair index. Only
 *  canonical structures (i.e. ones not containing lonely pairs) are generated.
 *  Only six special cases need to be checked. The resulting structures are
 *  pushed onto the global structure stack.
 */
void RNA2_move_noLP_bpshift() {
    /* Outside: shift i to left or j to right */
    int i, j, k;
    for( i=0; i<seq_len(); i++)
        if( is_opening_bp(i))        /* Handle each pair only once */
        {
            j = pairs_with(i);
            open_bp( i, j);
            if(is_bp(i+1, j-1) && !is_ins_lonely(i+1, j-1)) /* Outside & cross moves */
            {
                if( j < seq_len()-1)
                {
                    k = pairs_with(j+1);   /* i -- check only two possible pair (k+1,j)/(j,k+1) [cross] */
                    if( k >= 0
                        && k < seq_len()-1
                        && k != i-1
                        && !is_paired(k+1)
                        && valid_bp(k+1,j)
                    )
                    {
                        if( k<j)        /* Shift i to left */
                        {
                            close_bp( k+1, j);
                            if(VERB)
                                fprintf( stderr, "pushing sil %s\n", dot_bracket_str());
                        }
                        else            /* Cross i to right side */
                        {
                            close_bp( j, k+1);
                            if(VERB)
                                fprintf( stderr, "pushing sic %s j=%d k+1=%d\n", dot_bracket_str(), j, k+1);
                        }
                        push(dot_bracket_str());

                        if( k<j)
                            open_bp( k+1, j);
                        else
                            open_bp( j, k+1);
                    }
                }
                if( i > 0)
                {
                    k = pairs_with(i-1);  /* j -- check only two possible pair (i,k-1)/(k-1,i) [cross] */
                    if( k>0 && k!=j+1 && !is_paired(k-1) && valid_bp( i,k-1))
                    {
                        /*say "Shift j to right / cross: found candidate";*/
                        if( i<k)
                        {
                            close_bp( i, k-1);
                            if(VERB)
                                fprintf( stderr, "pushing sjr %s\n", dot_bracket_str());
                        }
                        else
                        {
                            close_bp( k-1, i);
                            if(VERB)
                                fprintf( stderr, "pushing sjc %s\n", dot_bracket_str());
                        }
                        push(dot_bracket_str());
                        if( i<k)
                            open_bp( i, k-1);
                        else
                            open_bp( k-1, i);
                    }
                }
            }
            /* Inside: shift i to right or j to left */
            if( i>0 && is_bp(i-1, j+1) && !is_out_lonely( i-1, j+1))
            {
                k = pairs_with(j-1);  /* i -- check only possible pair (k-1,j) */
                if( k>i+1 && !is_paired(k-1) && valid_bp( k-1, j)) /* Shift i to right */
                {
                    close_bp( k-1, j);
                    push(dot_bracket_str());
                    if(VERB)
                        fprintf( stderr, "pushing sir %s\n", dot_bracket_str());
                    open_bp( k-1, j);
                }
                k = pairs_with(i+1); /* j -- check only possible pair (i,k+1) */
                if( k>=0 && k<j-1 && !is_paired(k+1) && valid_bp(i,k+1))
                {
                    /* say "Shift j to left: found candidate"; */
                    close_bp( i, k+1);
                    push(dot_bracket_str());
                    if(VERB)
                        fprintf( stderr, "pushing sjl %s\n", dot_bracket_str());
                    open_bp( i, k+1);
                }
            }
            close_bp( i, j);
        }
}


/**
 * @brief Generate all neighbors by insertions of a base pair (i,j)
 *  where j>i, such that j is in the same loop as i (=>non-crossing)
 * @param i left position of base pairs
 * @param avoid_end right end, where no neighbor is generated
 * @note avoid_end is used when generating shift moves; for regular bp
 *  insertions, set avoid_end=i. The resulting structures are pushed onto the
 *  global structure stack.
 */
void RNA2_move_bpins_to_right(int i, int avoid_end)
{
    int j;
    for (j=i+1;j<seq_len(); j++) {

        /* jump over inner base pairs of current loop */
        while(j<seq_len() && is_opening_bp(j)) {
            j=pairs_with(j)+1;
        }

        if (is_closing_bp(j)) { /* here, the current loop ends */
            break;
        } else {
            /* j unpaired,
               register neighbor */
            if ( j-i>=MYTURN && j!=avoid_end && valid_bp(i,j)) {
                close_bp(i,j);
                if(VERB)
                    fprintf( stderr, "pushing %s %s\n", i==avoid_end ? "ins" : "shr", dot_bracket_str());
                push(dot_bracket_str());
                open_bp(i,j);
            }
        }
    }
}


/* analogous to RNA2_move_bpins_to_right */
void RNA2_move_bpins_to_left(int i, int avoid_end) {
    int j;
    for (j=i;j>0;) {
        j--;

        while(j>=0 && is_closing_bp(j)) {
            j=pairs_with(j);
            if (j==0)
                return;
            j--;
        }

        if (is_opening_bp(j)) {
            break;
        } else {
            /* j unpaired,
               register neighbor */
            if ( i-j>=MYTURN && j!=avoid_end && valid_bp(j,i)) {
                close_bp(j,i);
                if(VERB)
                    fprintf( stderr, "pushing %s %s\n",
                             i==avoid_end ? "ins" : "shl", dot_bracket_str());
                push(dot_bracket_str());
                open_bp(j,i);
            }
        }
    }
}


/**
 * @brief Generate neighbors by adding or deleting basepairs. Shift moves are
 *  not permitted, loneliness of pairs is ignored. The resulting structures are
 *  pushed onto the global structure stack.
 */
void RNA2_move_std() {
    unsigned int i;
    unsigned int j;

    /* generate base pair deletions */
    for (i=0; i<seq_len(); i++) {
        if (is_opening_bp(i)) {
            j = pairs_with(i);
            open_bp(i,j);
            if( VERB)
                fprintf( stderr, "pushing del %s\n", dot_bracket_str());
            push(dot_bracket_str());
            close_bp(i,j);
        }
    }

    /* generate base pair insertions */
    for (i=0; i<seq_len(); i++) {
        if (!is_paired(i))
            RNA2_move_bpins_to_right(i,i);
    }

    /* generate base pair shifts */
    if (r2d->shift) {
        int other_end;
        for (i=0; i<seq_len(); i++) {
            if (is_opening_bp(i) || is_closing_bp(i)) {
                other_end=pairs_with(i);

                open_bp(i,other_end);

                /* generate shifts to the right*/
                RNA2_move_bpins_to_right(i,other_end);

                /* generate shifts to the left */
                RNA2_move_bpins_to_left(i,other_end);

                /* restore basepair */
                if (other_end<i) {
                    close_bp(other_end,i);
                } else {
                    close_bp(i,other_end);
                }
            }
        }
    }
}


/**
 * @brief Parses the form string of the r2d structure and sets the pairtable
 *  accordingly.
 */
void parse_structure()
{
    int i;
    int siz=-1;

    for (i=0; i<seq_len(); i++) {
        if (is_opening_bp(i)) {
            r2d->stack[++siz]=i;
        } else if (is_closing_bp(i)) {
            if( !valid_bp( i, r2d->stack[siz]))
            {
                fprintf( stderr, "ERROR: invalid basepair (%c,%c) in"
                                 " parse_structure()\n",
                         r2d->seq[i], r2d->seq[r2d->stack[siz]]);
                exit(1);
            }
            r2d->pt[i]=r2d->stack[siz];
            r2d->pt[r2d->stack[siz]]=i;
            siz--;
        } else { /* no bracket => assume unpaired */
            r2d->pt[i] = UNPRD;
        }
    }
}


/**
 * @brief Generate all neighbors of the passed structure. The flags 'shift'
 * and 'noLP' of the r2d structure are taken into account.
 * @param structure Structure for which to generate neighbors
 */
void RNA2_move_it(char * structure) {
    if( strlen(structure) != seq_len() )
    {
        fprintf(stderr, "ERROR: length of sequence and passed structure "
                         "differ in RNA2_move_it()\n");
        exit(1);
    }
    strcpy(r2d->form, structure);
    parse_structure();  /* update pairtable */
    reset_stapel();     /* Clear global structure stack */

/*    if( VERB)
        fprintf( stderr, "                      1         2         3         4         5\n"
                     "            012345678901234567890123456789012345678901234567890\n"
                     "            %s shift=%d noLP=%d\n",
                 dot_bracket_str(), r2d->shift, r2d->noLP);
 */
    if (r2d->noLP) {
        RNA2_move_noLP_bpins();
        RNA2_move_noLP_bpdel();

        if (r2d->shift) {
            RNA2_move_noLP_bpshift();
        }
    } else {
        RNA2_move_std();
    }
    return;
}


void RNA2_free(void) {
    if(isInitialized)
    {
        isInitialized = 0;

        free(r2d->seq);
        free(r2d->form);
        free(r2d->pt);
        free(r2d->stack);
        free(r2d);
    }

    return;
}


char const* RNA2_get_seq()
{
    return r2d->seq;
}


const char *RNA2_get_form()
{
    return r2d->form;
}


int RNA2_get_len()
{
    return seq_len();
}

/*                                   EOF                                    */

