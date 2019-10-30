/***************************************************************************
          constraints.h  -  functions dealing with sequence constraints
                             -------------------
    begin                : 23-03-2006
    copyright            : (C) 2006 by Anke Busch
    email                : abusch@informatik.uni-freiburg.de
***************************************************************************/

#ifndef _CONSTRAINTS__
#define _CONSTRAINTS__

#include <stdlib.h>
#include "secis_basics.h"

using namespace std;

int Check_iu(char *secis_nuc);
int Check_secis_nuc_and_struct(char *secis_nuc, int *base_pairs);
int** getSeqConstraints(char *secis_nuc);

#endif   // _CONSTRAINTS_
