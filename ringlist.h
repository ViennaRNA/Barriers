/* Last changed Time-stamp: <2002-04-18 22:08:09 studla> */
/* ringlist.h */

#ifndef _ringlist_h
#define _ringlist_h

extern void RNA_init(char *sequence, int shift, int nolp);
extern void RNA_move_it(char *struc);
extern void RNA_free_rl(void);

extern void SPIN_move_it(char *sting);
extern void SPIN_complement_move_it(char *string);

extern void String_move_it(char *string);
extern void String_set_alpha(char *alpha);

extern void NNI_move_it(char *struc);

extern void Transpos_move_it(char *);
extern void CTranspos_move_it(char *);
extern void Reversal_move_it(char *);

extern void EXCH_move_it(char *);

extern char *pack_spin(const char *spin);
extern char *unpack_spin(const char *packed);

extern void LIST_move_it(char *);
extern void  put_ADJLIST(char *); 
#endif
