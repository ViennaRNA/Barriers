#define NBASES 8
static char Law_and_Order[]         = "_ACGUXKI";
static int  BP_pair[NBASES][NBASES] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

#define MAXALPHA 20       /* maximal length of alphabet */

static short  alias[MAXALPHA + 1];
static int    pair[MAXALPHA + 1][MAXALPHA + 1];
/* rtype[pair[i][j]]:=pair[j][i] */
static int    rtype[8] = {
  0, 2, 1, 4, 3, 6, 5, 7
};

#define ENCODE(C) ((strchr(Law_and_Order, \
                           (C)) == 0) ? 0 : (strchr(Law_and_Order, (C)) - Law_and_Order))

static void
make_pair_matrix(void)
{
  int i, j;

  for (i = 0; i < 5; i++)
    alias[i] = i;
  alias[5]  = 3;  /* X <-> G */
  alias[6]  = 2;  /* K <-> C */
  alias[7]  = 0;  /* I <-> default base '@' */
  for (i = 0; i < NBASES; i++)
    for (j = 0; j < NBASES; j++)
      pair[i][j] = BP_pair[i][j];
  /* if (noGU) pair[3][4] = pair[4][3] =0; */
}
