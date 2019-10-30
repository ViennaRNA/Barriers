#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "stapel.h"

#define NUM 5555

static char         *ALPHABET;
static int          ALPHASIZE;
static char         *ADJLIST;
/* static char *get_ADJLIST(void); */
static char **move;


void
initialize_crankshaft(void)
{
  int i;

  move = (char **)calloc(NUM, sizeof(char *));
  for (i = 0; i < NUM; i++)
    move[i] = (char *)calloc(5, sizeof(char));

  switch (ALPHASIZE) {
    /* 123 */
    case 3:                         /* FLR : SQ lattice */
      strcpy(move[232], "FLF\0");
      strcpy(move[121], "LRL\0");   /* LRL <-> FLF */
      strcpy(move[323], "FRF\0");
      strcpy(move[131], "RLR\0");   /* RLR <-> FRF */
      strcpy(move[231], "FLR\0");
      strcpy(move[123], "LRF\0");   /* LRF <-> FLR */
      strcpy(move[321], "FRL\0");
      strcpy(move[132], "RLF\0");   /* RLF <-> FRL */
      strcpy(move[2332], "RLLR\0");
      strcpy(move[3223], "LRRL\0"); /* LRRL <-> RLLR */
      break;
    case 5:                         /* 12345 */
                                    /* FLRUD: SC lattice */
      strcpy(move[231], "DDR\0");
      strcpy(move[553], "LRF\0");   /* LRF <-> DDR */
      strcpy(move[232], "DDD\0");
      strcpy(move[555], "LRL\0");   /* LRL <-> DDD */
      strcpy(move[233], "DDU\0");
      strcpy(move[554], "LRR\0");   /* LRR <-> DDU */
      strcpy(move[235], "DDF\0");
      strcpy(move[551], "LRD\0");   /* LRD <-> DDF */
      strcpy(move[321], "DUL\0");
      strcpy(move[542], "RLF\0");   /* RLF <-> DUL */
      strcpy(move[322], "DUU\0");
      strcpy(move[544], "RLL\0");   /* RLL <-> DUU */
      strcpy(move[323], "DUD\0");
      strcpy(move[545], "RLR\0");   /* RLR <-> DUD */
      strcpy(move[324], "DUF\0");
      strcpy(move[541], "RLU\0");   /* RLU <-> DUF */
      strcpy(move[451], "RRD\0");
      strcpy(move[335], "UDF\0");   /* UDF <-> RRD */
      strcpy(move[453], "RRF\0");
      strcpy(move[331], "UDR\0");   /* UDR <-> RRF */
      strcpy(move[454], "RRR\0");
      strcpy(move[333], "UDU\0");   /* UDU <-> RRR */
      strcpy(move[455], "RRL\0");
      strcpy(move[332], "UDD\0");   /* UDD <-> RRL */
      strcpy(move[441], "LLU\0");
      strcpy(move[224], "UUF\0");   /* UUF <-> LLU */
      strcpy(move[442], "LLF\0");
      strcpy(move[221], "UUL\0");   /* UUL <-> LLF */
      strcpy(move[444], "LLL\0");
      strcpy(move[222], "UUU\0");   /* UUU <-> LLL */
      strcpy(move[445], "LLR\0");
      strcpy(move[223], "UUD\0");   /* UUD <-> LLR */
      strcpy(move[121], "LUL\0");
      strcpy(move[242], "FLF\0");   /* FLF <-> LUL */
      strcpy(move[122], "LUU\0");
      strcpy(move[244], "FLL\0");   /* FLL <-> LUU */
      strcpy(move[123], "LUD\0");
      strcpy(move[245], "FLR\0");   /* FLR <-> LUD */
      strcpy(move[124], "LUF\0");
      strcpy(move[241], "FLU\0");   /* FLU <-> LUF */
      strcpy(move[421], "FUL\0");
      strcpy(move[142], "ULF\0");   /* ULF <-> FUL */
      strcpy(move[422], "FUU\0");
      strcpy(move[144], "ULL\0");   /* ULL <-> FUU */
      strcpy(move[423], "FUD\0");
      strcpy(move[145], "ULR\0");   /* ULR <-> FUD */
      strcpy(move[424], "FUF\0");
      strcpy(move[141], "ULU\0");   /* ULU <-> FUF */
      strcpy(move[151], "DRD\0");
      strcpy(move[535], "FDF\0");   /* FDF <-> DRD */
      strcpy(move[153], "DRF\0");
      strcpy(move[531], "FDR\0");   /* FDR <-> DRF */
      strcpy(move[154], "DRR\0");
      strcpy(move[533], "FDU\0");   /* FDU <-> DRR */
      strcpy(move[155], "DRL\0");
      strcpy(move[532], "FDD\0");   /* FDD <-> DRL */
      strcpy(move[131], "RDR\0");
      strcpy(move[353], "FRF\0");   /* FRF <-> RDR */
      strcpy(move[132], "RDD\0");
      strcpy(move[355], "FRL\0");   /* FRL <-> RDD */
      strcpy(move[133], "RDU\0");
      strcpy(move[354], "FRR\0");   /* FRR <-> RDU */
      strcpy(move[135], "RDF\0");
      strcpy(move[351], "FRD\0");   /* FRD <-> RDF */
      strcpy(move[3524], "LURD\0");
      strcpy(move[2435], "RDLU\0"); /* RDLU <-> LURD */
      strcpy(move[5342], "ULDR\0");
      strcpy(move[4253], "DRUL\0"); /* DRUL <-> ULDR */
      break;
    default:
      break;
  }
  return;
}


void
Q_mem_cleanup(void)
{
  int i;

  for (i = 0; i < NUM; i++)
    free(move[i]);
  free(move);
}


void
String_move_it_crankshaft(char *string)
{
  int   i, j, length, rep_len, start = 1, found_em = 0;
  char  cc, *s, *tr_str;

  s = strdup(string);

  length  = strlen(string);
  rep_len = 4; /* permute max. 4 letters for crankshaft moves */
  tr_str  = (char *)calloc(rep_len + 1, sizeof(char));

  while (start <= (length - rep_len) + 1) {
    int zahl;
    for (i = 0; i < rep_len; i++) {
      switch (s[start + i]) {
        case 'F':
          strcat(tr_str + i, "1");
          break;
        case 'L':
          strcat(tr_str + i, "2");
          break;
        case 'R':
          strcat(tr_str + i, "3");
          break;
        case 'U':
          strcat(tr_str + i, "4");
          break;
        case 'D':
          strcat(tr_str + i, "5");
          break;
        default:
          break;
      }
      if (i >= 2) {
        /* > 3 letters translated */
        zahl = atoi(tr_str);
        if (strpbrk(move[zahl], ALPHABET)) {
          for (j = 0; j <= i; j++)
            s[start + j] = move[zahl][j];
          /*  fprintf(stdout, "SS %s : pushed %s onto stack\n", string, s); */
          push(s);
          free(s);
          s = strdup(string);
          memset(tr_str, 0, rep_len + 1);
          i++;  /* found a 3-letter replacement, no need */
                /* to find a 4-letter replacement too */
        }
      }
    }
    memset(tr_str, 0, rep_len + 1);
    start++;
  }

  free(s);
  if (tr_str != NULL)
    free(tr_str);

  /* second: manipulate end of the string */
  s = strdup(string);
  switch (ALPHASIZE) {
    case 3: /* SQ */
      if (strcmp(string + (length - 2), "LR") == 0) {
        strcpy(s + (length - 2), "FL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "FL") == 0) {
        strcpy(s + (length - 2), "LR\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "FR") == 0) {
        strcpy(s + (length - 2), "RL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "RL") == 0) {
        strcpy(s + (length - 2), "FR\0");
        found_em = 1;
      }

      break;
    case 5: /* SC */
      if (strcmp(string + (length - 2), "FL") == 0) {
        strcpy(s + (length - 2), "LU\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "LU") == 0) {
        strcpy(s + (length - 2), "FL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "FR") == 0) {
        strcpy(s + (length - 2), "RD\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "RD") == 0) {
        strcpy(s + (length - 2), "FR\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "FU") == 0) {
        strcpy(s + (length - 2), "UL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "UL") == 0) {
        strcpy(s + (length - 2), "FU\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "FD") == 0) {
        strcpy(s + (length - 2), "DR\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "DR") == 0) {
        strcpy(s + (length - 2), "FD\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "LL") == 0) {
        strcpy(s + (length - 2), "UU\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "UU") == 0) {
        strcpy(s + (length - 2), "LL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "LR") == 0) {
        strcpy(s + (length - 2), "DD\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "DD") == 0) {
        strcpy(s + (length - 2), "LR\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "RL") == 0) {
        strcpy(s + (length - 2), "DU\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "DU") == 0) {
        strcpy(s + (length - 2), "RL\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "RR") == 0) {
        strcpy(s + (length - 2), "UD\0");
        found_em = 1;
      } else if (strcmp(string + (length - 2), "UD") == 0) {
        strcpy(s + (length - 2), "RR\0");
        found_em = 1;
      }

      break;
    default:
      break;
  }

  if (found_em)
    push(s);

  /*   fprintf(stdout, "SS %s : pushed %s onto stack (end)\n", string, s); */

  free(s);

  /* third: do end move(s): derived from pivot-routine */
  s   = strdup(string);
  cc  = s[length - 1]; /* the last relative move */
  for (j = 0; j < ALPHASIZE; j++) {
    /* exchange current char with all possible ones */
    if (cc == ALPHABET[j])
      continue;                     /* skip current char, needs not to be changed */

    s[length - 1] = ALPHABET[j];
    push(s);
    /*  fprintf(stdout, "SS %s : pushed %s onto stack (end pivot)\n", string, s); */
  }
  free(s);
  return;
}


void
String_move_it(char *string)
{
  int   i, j, length, ID;
  char  *s;

  length  = strlen(string);
  s       = (char *)space((length + 1) * sizeof(char));

  for (i = 0; i < length; i++) {
    ID = 0;
    for (j = 0; j < ALPHASIZE; j++) {
      if (ALPHABET[j] != string[i]) {
        strcpy(s, string);
        s[i] = ALPHABET[j];
        push(s);
      } else {
        ID++;
      }
    }
    if (ID != 1)
      fprintf(stderr,
              "Warning: Input string not from given Alphabet\n");
  }
  free(s);
}


void
String_set_alpha(char *alpha)
{
  ALPHABET  = alpha;
  ALPHASIZE = strlen(alpha);
}


void
SPIN_move_it(char *string)
{
  /* generate 1-point error mutants */
  int   i, length;
  char  *s;

  length  = strlen(string);
  s       = (char *)space((length + 1) * sizeof(char));

  for (i = 0; i < length; i++) {
    strcpy(s, string);
    if (s[i] == '+')
      s[i] = '-';
    else
      s[i] = '+';

    push(s);
  }
  free(s);
}


void
SPIN_complement_move_it(char *string)
{
  /* complement string after position k for all k */
  int   i, length;
  char  *s;

  length  = strlen(string);
  s       = strdup(string);
  for (i = length - 1; i >= 0; i--) {
    if (s[i] == '+')
      s[i] = '-';
    else
      s[i] = '+';

    push(s);
  }
  free(s);
}


/* Move-sets on Permutations */
static int *
String2Perm(char *perm)
{
  int   i, k, n;
  int   *P;
  char  *s, *line;

  line = strdup(perm);
  for (k = i = 0; i < strlen(line); i++)
    if (line[i] == ',')
      k++;

  n     = k + 1;
  P     = (int *)space(sizeof(int) * (n + 1));
  P[0]  = n;
  s     = strtok(line, ",\n");
  for (i = 1; i <= n; i++) {
    (void)sscanf(s, "%d", &k);
    P[i]  = k;
    s     = strtok(NULL, ",\n");
  }
  free(line);
  return P;
}


static char *
Perm2String(int *P)
{
  char  tmp[1000];
  char  digit[10];
  char  *string;
  int   i, n;

  n       = P[0];
  tmp[0]  = '\0';
  for (i = 1; i < n; i++) {
    sprintf(digit, "%d,", P[i]);
    strcat(tmp, digit);
  }
  sprintf(digit, "%d", P[i]);
  strcat(tmp, digit);
  string = strdup(tmp);
  return string;
}


void
Transpos_move_it(char *perm)
{
  int   i, j, jj, n;
  int   *P;
  char  *s;

  P = String2Perm(perm);
  n = P[0];

  for (i = 1; i < n; i++)
    for (j = i + 1; j <= n; j++) {
      jj    = P[i];
      P[i]  = P[j];
      P[j]  = jj;                   /* transpose i,j */
      s     = Perm2String(P);
      jj    = P[i];
      P[i]  = P[j];
      P[j]  = jj;                   /* transpose i,j again to restore P */
      push(s);
      free(s);
    }
  free(P);
}


void
CTranspos_move_it(char *perm)
{
  int   i, j, jj, n;
  int   *P;
  char  *s;

  P = String2Perm(perm);
  n = P[0];

  for (i = 1; i < n; i++) {
    j     = i + 1;
    jj    = P[i];
    P[i]  = P[j];
    P[j]  = jj;                   /* transpose i,j */
    s     = Perm2String(P);
    jj    = P[i];
    P[i]  = P[j];
    P[j]  = jj;                   /* transpose i,j again to restore P */
    push(s);
    free(s);
  }
  free(P);
}


void
Reversal_move_it(char *perm)
{
  int   i, j, k, n;
  int   *P, *P1;
  char  *s;

  P   = String2Perm(perm);
  n   = P[0];
  P1  = space(sizeof(int) * (n + 1));

  for (i = 1; i < n; i++)
    for (j = i + 1; j <= n; j++) {
      for (k = 1; k < i; k++)
        P1[k] = P[k];
      for (k = 0; k < j - i; k++)
        P1[i + k] = P[j - k];
      for (k = j + 1; k <= n; k++)
        P1[k] = P[k];
      s = Perm2String(P);
      push(s);
      free(s);
    }
  free(P);
  free(P1);
}


/* Move-sets for Trees */
#include "trees.h"

void
NNI_move_it(char *string)
{
  /* Nearest neighbor Interchange moves */
  int   i, nl, nneigh;
  char  *s;
  Tree  *T;
  Tree  *T_NNI;

  /* printf("%s\n",string); */

  nl      = number_of_leaves(string);
  nneigh  = 2 * (nl - 3);
  T       = string2Tree(string);
  T_NNI   = Make_all_NNI(T[0]);

  for (i = 0; i < nneigh; i++) {
    s = Tree2string(T_NNI[i]);
    push(s);
    free(s);
  }

  FreeTree(T, 1);
  FreeTree(T_NNI, nneigh);
}


static int spin_len;

char *unpack_spin(const char *packed);


char *
pack_spin(const char *spin)
{
  int           mask[7] = {
    64, 32, 16, 8, 4, 2, 1
  };
  int   i, j, k, l;
  char  *packed;

  spin_len  = strlen(spin);
  l         = (spin_len + 6) / 7;
  packed    = (char *)space(l * sizeof(char) + 1);
  for (i = j = 0; i < spin_len; j++) {
    for (k = 0; (k < 7) && (i < spin_len); k++, i++) {
      if (spin[i] == '+')
        packed[j] |= mask[k];
      else if (spin[i] != '-')
        fprintf(stderr, "Junk in spin %s\n", spin);
    }
    packed[j]++;
  }
#if 0
  {
    char *s;
    s = unpack_spin(packed);
    if (strcmp(s, spin) != 0)
      fprintf(stderr, "Error in pack_spin %s %s %s\n", spin, s, packed);

    free(s);
  }
#endif
  return packed;
}


char *
unpack_spin(const char *packed)
{
  int   i, j, k, l;
  int   mask[7] = {
    64, 32, 16, 8, 4, 2, 1
  };
  char  *spin;

  l     = strlen(packed);
  spin  = space((7 * l + 1) * sizeof(char));
  for (i = j = 0; j < l; j++) {
    int p;
    p = packed[j] - 1;
    for (k = 0; k < 7; k++)
      spin[i++] = (p & mask[k]) ? '+' : '-';
  }
  spin[spin_len] = '\0';
  return spin;
}


void
EXCH_move_it(char *string)
{
  int   i, j, length;
  char  *s;

  length  = strlen(string);
  s       = (char *)space((length + 1) * sizeof(char));

  for (i = 0; i < length; i++) {
    for (j = 0; j < length; j++) {
      if ((string[i] == '+') && (string[j] == '-')) {
        strcpy(s, string);
        s[i]  = '-';
        s[j]  = '+';
        push(s);
      }
    }
  }
  free(s);
}


/********************************************************************/

void
put_ADJLIST(char *A)
{
  if (ADJLIST != NULL)
    free(ADJLIST);

  ADJLIST = strdup(A);
}


void
LIST_move_it(char *string)
{
  char *s, *token;

  /* parse adjacency list and push the token on the stack */
  /* using strtok */

  token = strtok(ADJLIST, ":");
  while (token) {
    s = strdup(token);
    push(s);
    token = strtok(NULL, ":");
  }
  free(ADJLIST);
  ADJLIST = NULL;
}


/********************************************************************/
