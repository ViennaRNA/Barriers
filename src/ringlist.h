/* ringlist.h */

#ifndef BARRIERS_RINGLIST_H
#define BARRIERS_RINGLIST_H

extern int MYTURN;

void
RNA_init(char *sequence,
         int  shift,
         int  nolp);


void
RNA_move_it(char *struc);


void
RNA_move_it_rates(char *form);


void
RNA_free_rl(void);


void
SPIN_move_it(char *sting);


void
SPIN_complement_move_it(char *string);


void
String_move_it(char *string);


void
String_move_it_crankshaft(char *string);


void
String_set_alpha(char *alpha);


void
initialize_crankshaft(void);


void
Q_mem_cleanup(void);


void
NNI_move_it(char *struc);


void
Transpos_move_it(char *);


void
CTranspos_move_it(char *);


void
Reversal_move_it(char *);


void
EXCH_move_it(char *);


char *
pack_spin(const char *spin);


char *
unpack_spin(const char *packed);


void
LIST_move_it(char *);


void
put_ADJLIST(char *);


#endif
