/* Last changed Time-stamp: <2001-06-11 14:09:40 studla> */
/* ringlist.h */

#ifndef _ringlist_h
#define _ringlist_h

extern void RNA_init(char *sequence, int shift, int nolp);
extern void RNA_move_it(char *struc);
extern void RNA_free_rl(void);

extern void SPIN_move_it(char *struc);
extern void String_move_it(char *string);
extern void String_set_alpha(char *alpha);

extern void NNI_move_it(char *struc);

extern void Transpos_move_it(char *);
extern void CTranspos_move_it(char *);
extern void Reversal_move_it(char *);

extern void EXCH_move_it(char *);

extern char *pack_spin(const char *spin);
extern char *unpack_spin(const char *packed);
#endif
