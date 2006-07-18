/*
  Last changed Time-stamp: <2006-07-13 11:37:49 xtof>
  $Id: secis_neighbors.h,v 1.1 2006/07/18 14:08:49 xtof Exp $
*/

#ifndef _SECIS_NEIGHBORS_H_
#define _SECIS_NEIGHBORS_H_

#if __cplusplus
#  define BEGIN_C_DECLS extern "C" {
#  define END_C_DECLS   }
#else
#  define BEGIN_C_DECLS
#  define END_C_DECLS
#endif

BEGIN_C_DECLS
void initialize_SECIS(char* s_nuc, char* s_struct, char* p_seq,
		      int max_m, int min_as);

void SECIS_move_it(char* rna_seq);
END_C_DECLS

#endif // _SECIS_NEIGHBORS_H_

/* End of file */

