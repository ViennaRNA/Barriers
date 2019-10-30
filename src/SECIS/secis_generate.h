/***************************************************************************
                generate.h  -  generating valid sequences
                             -------------------
    begin                : 27-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/

#ifndef _GENERATE__
#define _GENERATE__

#include <stdlib.h>
#include "secis_basics.h"

using namespace std;

int* init_recursion(char* s_nuc, char* s_struct, char* p_seq,
		    int max_m, int min_as);
void gen_seq(int pos, int* seq);


#endif   // _GENERATE_

