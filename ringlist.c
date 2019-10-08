/* Last changed Time-stamp: <2017-10-30 14:22:26 mtw> */
/* ringlist.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "utils.h"
#include "pair_mat.h"
#include "stapel.h"

static char UNUSED  rcsid[] = "$Id: ringlist.c,v 1.1 2001/04/05 08:00:57 ivo Exp $";

int                 MYTURN = 4;

typedef struct _rlItem {
  int             nummer; /* number of base in sequence */
  char            typ;    /* basepair q or not p - r virtual root */
  unsigned short  base;   /* 0<->unknown, 1<->A, 2<->C, 3<->G, 4<->U */
  struct _rlItem  *up;
  struct _rlItem  *next;
  struct _rlItem  *prev;
  struct _rlItem  *down;
}rlItem;

static char   *form     = NULL; /* array 4 (.)-structure */
static char   *farbe    = NULL; /* array 4 sequence */
static int    len       = 0;    /* length of sequence */
static int    poListop  = 0;    /* polist counter = no_of_bp */
static rlItem *rl       = NULL; /* array 4 ringlist */
static rlItem *wurzl    = NULL; /* virtual root of ringlist-tree */
static rlItem **poList  = NULL; /* post order list of bp's */

static int    xtof;             /* do shift moves */
static int    noLP;             /* no lonely pairs move-set */
static void ini_or_reset_rl(char  *seq,
                            char  *struc);


/* public functiones */
void RNA_init(char  *seq,
              int   xtof,
              int   noLP);


void RNA_move_it(char *struc);


void RNA_move_it_rates(char *form);


void RNA_free_rl(void);


#ifdef HARDCORE_DEBUG
void rl_status(void);


#endif

/* private functiones */
static void struc2tree(char *struc);


static void open_bp(rlItem *i);


static void close_bp(rlItem *i,
                     rlItem *j);


static void inb(rlItem *root);


static void inb_nolp(rlItem *root);


static void dnb(rlItem *rli);


static void dnb_nolp(rlItem *rli);


static void fnb(rlItem *rli);


static void make_poList(rlItem *root);


void
RNA_init(char *seq,
         int  shift,
         int  nolp)
{
  xtof  = shift;
  noLP  = nolp;
  farbe = seq;
}


/**/
static void
ini_or_reset_rl(char  *seq,
                char  *struc)
{
  int   i;
  char  *pos;

  len = strlen(seq);
  if (wurzl == NULL) {
    form  = strdup(struc);
    form  = (char *)xrealloc(form, (strlen(struc) + 3) * sizeof(char));
    /* form[len+1] = '\0'; */
    form[strlen(struc)]     = '\0';
    form[strlen(struc) + 1] = '\0';
    form[strlen(struc) + 2] = '\0';
    farbe                   = strdup(seq);
    /*      update_fold_params(); */
    make_pair_matrix();
    poList              = (rlItem **)calloc(len, sizeof(rlItem *));
    rl                  = (rlItem *)calloc(len + 1, sizeof(rlItem));
    wurzl               = (rlItem *)calloc(1, sizeof(rlItem));
    wurzl->typ          = 'r';
    wurzl->nummer       = -1;
    wurzl->down         = &rl[len];
    poList[poListop++]  = wurzl;

    for (i = 0; i < len; i++) {
      rl[i].typ = 'u';
      pos       = strchr(Law_and_Order, seq[i]);
      if (pos == NULL)
        rl[i].base = 0;
      else
        rl[i].base = pos - Law_and_Order;

      rl[i].nummer  = i;
      rl[i].next    = &rl[i + 1];
      rl[i].prev    = ((i == 0) ? &rl[len] : &rl[i - 1]);
      rl[i].up      = rl[i].down = NULL;
    }
    rl[i].nummer  = i;
    rl[i].base    = 0;
    rl[i].next    = &rl[0];     /* rl.next ist jetzt kreis */
    rl[i].prev    = &rl[i - 1]; /* rl.prev ist jetzt kreis */
    rl[i].up      = wurzl;
    rl[i].typ     = 'x';
    /* ini_stapel(len); */
  } else {
    /* reset ringlist */
    strcpy(form, struc);
    for (i = 0; i < len; i++) {
      rl[i].typ   = 'u';
      rl[i].next  = &rl[i + 1];
      rl[i].prev  = ((i == 0) ? &rl[len] : &rl[i - 1]);
      rl[i].up    = rl[i].down = NULL;
    }
    rl[i].next  = &rl[0];     /* rl.next ist jetzt kreis */
    rl[i].prev  = &rl[i - 1]; /* rl.prev ist jetzt kreis */
    rl[i].up    = wurzl;
  }

  struc2tree(struc);
  reset_stapel();
}


/**/
void
RNA_free_rl(void)
{
  free(rl);
  free(wurzl);
  free(farbe);
  free(form);
  free(poList);
}


/**/
static void
struc2tree(char *struc)
{
  int     ipos, jpos, balance = 0;
  char    *struc_copy;
  rlItem  *rli, *rlj;

  struc_copy = strdup(struc);
  for (ipos = 0; ipos < len; ipos++) {
    if (struc_copy[ipos] == ')') {
      jpos              = ipos;
      struc_copy[ipos]  = '.';
      balance++;
      while (struc_copy[--ipos] != '(');
      struc_copy[ipos] = '.';
      balance--;
      rli = &rl[ipos];
      rlj = &rl[jpos];
      close_bp(rli, rlj);
    }
  }
  if (balance) {
    fprintf(stderr, "startS2tree(): structure is not balanced !\n%s\n", struc);
    exit(1);
  }

  poListop = 1;
  make_poList(NULL);
  free(struc_copy);
}


/* close a particular base pair */
void
close_bp(rlItem *i,
         rlItem *j)
{
  rlItem *jn; /* points to j->next */

  /* change string representation */
  form[i->nummer] = '(';
  form[j->nummer] = ')';

  jn            = j->next;
  i->typ        = 'p';
  j->typ        = 'q';
  i->down       = j;
  j->up         = i;
  i->next->prev = j;
  j->next->prev = i;
  j->next       = i->next;
  i->next       = jn;
}


/* open a particular base pair */
void
open_bp(rlItem *i)
{
  rlItem *in; /* points to i->next */

  /* change string representation */
  form[i->nummer]       = '.';
  form[i->down->nummer] = '.';

  in            = i->next;
  i->typ        = 'u';
  i->down->typ  = 'u';
  i->next       = i->down->next;
  i->next->prev = i;
  in->prev      = i->down;
  i->down->next = in;
  i->down       = in->prev->up = NULL;
}


/* for a given tree, generate postorder-list */
static void
make_poList(rlItem *root)
{
  rlItem *stop, *rli;

  if (!root)
    root = wurzl;

  stop = root->down;
  for (rli = stop->next; rli != stop; rli = rli->next) {
    if (rli->typ == 'p') {
      /* wenn basenpaar untersuche subrlItem */
      poList[poListop++] = rli;
      make_poList(rli);
    }
  }
  return;
}


/* for a given tree, generate all neighbours according to the moveset */
void
RNA_move_it(char *form)
{
  int i;

  /* ini_or_reset_rl() allocates (length(seq)+2) characters space for
   * the structure, hence the ringlist-based neighbor routines tolerate
   * one additional char at the end of the structure, as e.g. in ligand
   * case */
  ini_or_reset_rl(farbe, form);
#ifdef DEBUG_NB
  fprintf(stderr, "m %s\n", form);
#endif
  if (noLP) {
    /* canonic neighbours only */
    for (i = 0; i < poListop; i++) {
      inb_nolp(poList[i]);
      if (i > 0)     /* virtual root should never be deleted or fliped */
        dnb_nolp(poList[i]);
        /*  if(xtof) fnb(poList[i]); */
    }
  } else {
    /* all neighbours */
    for (i = 0; i < poListop; i++) {
      inb(poList[i]);
      if (i > 0) {
        /* virtual root should never be deleted or fliped */
        dnb(poList[i]);
        if (xtof)
          fnb(poList[i]);
      }
    }
  }
}


/* special version of RNA_move_it that does additional moves from
 * starred to unstarred and vice versa */
void
RNA_move_it_rates(char *form)
{
  int   i, formlen;
  bool  hasstar = false;

#ifdef DEBUG_NB
  fprintf(stderr, "mR%s\n", form);
#endif
  formlen = strlen(form);
  if (form[formlen - 1] == '*')
    hasstar = true;

  RNA_move_it(form);

  /* now add special moves: starred -> unstarred and vice versa */
  if (hasstar == true) {
    /* for starred structure add unstarred neighbor */
    form[formlen - 1] = '\0';
    push(form);
    form[formlen - 1] = '*';
    form[formlen]     = '\0';
  } else {
    /* for an un-starred structure add a starred neighbor */
    form[formlen]     = '*';
    form[formlen + 1] = '\0';
    push(form);
    form[formlen] = '\0';
  }
}


/* for a given ringlist, generate all insert moves */
static void
inb(rlItem *root)
{
  rlItem *stop, *rli, *rlj;

  stop = root->down;
  for (rli = stop->next; rli != stop; rli = rli->next) {
    if (rli->typ == 'p')
      continue;

    for (rlj = rli->next; rlj != stop; rlj = rlj->next) {
      if (rlj->nummer - rli->nummer < MYTURN)
        continue;

      if (rlj->typ == 'p')
        continue;

      if (pair[rli->base][rlj->base]) {
        close_bp(rli, rlj);
#ifdef DEBUG_NB
        fprintf(stderr, "i%s\n", form);
#endif
        push(form);
        open_bp(rli);
      }
    }
  }
}


static void
inb_nolp(rlItem *root)
{
  rlItem *stop, *rli, *rlj;

  stop = root->down;
  for (rli = stop->next; rli != stop; rli = rli->next) {
    if (rli->typ == 'p')
      continue;

    for (rlj = rli->next; rlj != stop; rlj = rlj->next) {
      if (rlj->nummer - rli->nummer < MYTURN)
        continue;

      if (rlj->typ == 'p')
        continue;

      if (pair[rli->base][rlj->base]) {
        if (((rli->prev == stop && rlj->next == stop) && stop->typ != 'x') ||
            (rli->next == rlj->prev)) {
          /* base pair extends helix */
          close_bp(rli, rlj);
#ifdef DEBUG_NB
          fprintf(stderr, "iS%s\n", form);
#endif
          push(form);
          open_bp(rli);
        } else if ((rlj->nummer - rli->nummer >= MYTURN + 2) &&
                   (rli->next->typ != 'p' && rlj->prev->typ != 'p') &&
                   (rli->next->next != rlj->prev->prev) &&
                   (pair[rli->next->base][rlj->prev->base])) {
          /* double insert */
          close_bp(rli->next, rlj->prev);
          close_bp(rli, rlj);
          if (form[len] == '*') {
            form[len + 1] = 'D';
            form[len + 2] = '\0';
          } else {
            form[len]     = 'D';
            form[len + 1] = '\0';
          }

#ifdef DEBUG_NB
          fprintf(stderr, "iD%s\n", form);
#endif
          push(form);
          if (form[len] == '*') /* reset inserted D */
            form[len + 1] = '\0';
          else
            form[len] = '\0';

          /* form[len]='\0'; */
          open_bp(rli);
          open_bp(rli->next);
        }
      }
    }
  }
}


#if 0
/* for a given ringlist, generate all canonic inserte moves */
static void
inb_noLP(rlItem *root)
{
  rlItem *stop, *rli, *rlj;

  stop = root->down; /* j-pos of bp */
  /* loop ringlist 4 potential i-pos of bp */
  for (rli = rlj = stop->next; rli != stop; rli = rli->next) {
    /* helix elongation */
    if (rli->typ == 'p') {
      if (rli->prev == stop || rli->next == stop)
        continue;

      if (rli->prev->typ == 'p' || rli->next->typ == 'p')
        continue;

      if (pair[rli->prev->base][rli->next->base]) {
        /* enongate helix 2 the exterior */
        close_bp(rli->prev, rli->next);
        push(form);
        open_bp(rli->prev->up);
      }

      if (rli->down->next == rli->down || rli->down->prev == rli->down)
        continue;

      if (rli->down->next->typ == 'p' || rli->down->prev->typ == 'p')
        continue;

      if (pair[rli->down->next->base][rli->down->prev->base]) {
        /* enongate helix 2 the interior */
        close_bp(rli->down->next, rli->down->prev);
        push(form);
        open_bp(rli->down->next);
      }

      continue;
    }

    /* loop ringlist 4 potential j-pos of bp */
    for (rlj = rli->next; rlj != stop; rlj = rlj->next) {
      /* continue if bp violates minloop-condition */
      if (rlj->nummer - rli->nummer < MYTURN)
        continue;

      /* continue if j-pos is paired */
      if (rlj->typ == 'p')
        continue;

      /* found a legal bp */
      if (pair[rli->base][rlj->base]) {
        /* continue if immediate interior bp violates minloop-condition */
        if (rlj->prev->nummer - rli->next->nummer < MYTURN)
          continue;

        /* cont if i or j pos of immediate bp is allready paired */
        if (rli->next->typ == 'p' || rlj->prev->typ == 'p')
          continue;

        /* found immediate interior bp */
        if (pair[rli->next->base][rlj->prev->base]) {
          close_bp(rli->next, rlj->prev);
          close_bp(rli, rlj);
          push(form);
          open_bp(rli);
          open_bp(rli->next);
        }
      }
    }
  }
}


#endif

/* for a given ringlist, generate all delete moves */
static void
dnb(rlItem *rli)
{
  rlItem *rlj;

  rlj = rli->down;
  open_bp(rli);
  /* fprintf(stderr, "d%s\n", form); */
  push(form);
  close_bp(rli, rlj);
}


/* for a given ringlist, generate all canonic delete moves */
static void
dnb_nolp(rlItem *rli)
{
  rlItem  *rlj  = NULL;
  rlItem  *rlin = NULL; /* pointers to following pair in helix, if any */
  rlItem  *rljn = NULL;
  rlItem  *rlip = NULL; /* pointers to preceding pair in helix, if any */
  rlItem  *rljp = NULL;

  rlj = rli->down;
  if (rlj->next == rlj->prev) {
    /* immediate interior bp ? */
    rlin  = rlj->next;
    rljn  = rlin->down;
  }

  if (rli->prev == rli->next && rli->next->typ != 'x') {
    /* immediate exterior */
    rlip  = rli->next->up;
    rljp  = rli->next;
  }

  if (rlip == NULL && rlin && rljn->next != rljn->prev) {
    /* doubledelete */
    open_bp(rli);
    open_bp(rlin);
    if (form[len] == '*') {
      form[len + 1] = 'D';
      form[len + 2] = '\0';
    } else {
      form[len]     = 'D';
      form[len + 1] = '\0';
    }

#ifdef DEBUG_NB
    fprintf(stderr, "dD%s\n", form);
#endif
    push(form);
    if (form[len] == '*') /* reset inserted D */
      form[len + 1] = '\0';
    else
      form[len] = '\0';

    /* form[len]='\0'; */
    close_bp(rlin, rljn);
    close_bp(rli, rlj);
  } else {
    /* the following will work only if boolean expr are shortcicuited */
    if (rlip == NULL || (rlip->prev == rlip->next && rlip->prev->typ != 'x')) {
      if (rlin == NULL || (rljn->next == rljn->prev)) {
        open_bp(rli);
#ifdef DEBUG_NB
        fprintf(stderr, "dS%s\n", form);
#endif
        push(form);
        close_bp(rli, rlj);
      }
    }
  }
}


/* for a given ringlist, generate all flip moves */
static void
fnb(rlItem *rli)
{
  int     x;
  rlItem  *rlj, *stop, *help_rli, *help_rlj;

  stop = rli->down;
  /* examin interior loop of bp(ij) */
  for (rlj = stop->next; rlj != stop; rlj = rlj->next) {
    if ((rlj->typ == 'p') || (rlj->typ == 'q'))
      continue;                                    /* prevent shifting to paired position */

    if ((rlj->nummer - rli->nummer >= MYTURN) && (pair[rli->base][rlj->base])) {
      /* (ij)->(ik) i<k<j */
      open_bp(rli);         /* open original basepair */
      close_bp(rli, rlj);   /* close shifted basepair */
      push(form);
      open_bp(rli);         /* open shifted basepair */
      close_bp(rli, stop);  /* restore original basepair */
    }

    if ((stop->nummer - rlj->nummer >= MYTURN) && (pair[stop->base][rlj->base])) {
      /* (ij)->(kj) i<k<j */
      open_bp(rli);         /* open original basepair */
      close_bp(rlj, stop);  /* close shifted basepair */
      push(form);
      open_bp(rlj);         /* open shifted basepair */
      close_bp(rli, stop);  /* restore original basepair */
    }
  }
  /* examin exterior loop of bp(ij) */
  for (rlj = rli->next; rlj != rli; rlj = rlj->next) {
    if ((rlj->typ == 'p') || (rlj->typ == 'q') || (rlj->typ == 'x'))
      continue;

    x = rlj->nummer - rli->nummer;
    if (x < 0)
      x = -x;

    if ((x >= MYTURN) && (pair[rli->base][rlj->base])) {
      if (rli->nummer < rlj->nummer) {
        help_rli  = rli;
        help_rlj  = rlj;
      } else {
        help_rli  = rlj;
        help_rlj  = rli;
      }

      open_bp(rli);                 /* open original basepair */
      close_bp(help_rli, help_rlj); /* close shifted basepair */
      push(form);
      open_bp(help_rli);            /* open shifted basepair */
      close_bp(rli, stop);          /* restore original basepair */
    }

    x = rlj->nummer - stop->nummer;
    if (x < 0)
      x = -x;

    if ((x >= MYTURN) && (pair[stop->base][rlj->base])) {
      if (stop->nummer < rlj->nummer) {
        help_rli  = stop;
        help_rlj  = rlj;
      } else {
        help_rli  = rlj;
        help_rlj  = stop;
      }

      open_bp(rli);                 /* open original basepair */
      close_bp(help_rli, help_rlj); /* close shifted basepair */
      push(form);
      open_bp(help_rli);            /* open shifted basepair */
      close_bp(rli, stop);          /* restore original basepair */
    }
  }
}


#ifdef HARDCORE_DEBUG
/**/
void
rl_status(void)
{
  int i;

  printf("\n%s\n%s\n", farbe, form);
  for (i = 0; i <= len; i++) {
    printf("%2d %c %c %2d %2d %2d %2d\n",
           rl[i].nummer,
           i == len ? 'X' : farbe[i],
           rl[i].typ,
           rl[i].up == NULL ? 0 : (rl[i].up)->nummer,
           rl[i].down == NULL ? 0 : (rl[i].down)->nummer,
           (rl[i].prev)->nummer,
           (rl[i].next)->nummer);
  }
  printf("---\n");
}


#endif
