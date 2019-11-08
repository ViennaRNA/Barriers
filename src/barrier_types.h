#ifndef BARRIER_TYPES_H
#define BARRIER_TYPES_H

/* global structures */
#include "hash_table_linear_probing_lists/hash_tables.h"
typedef struct {
  unsigned long father;       /* which lmin do I merge with */
  char          *saddle;      /* structure of saddle point */
  float         E_saddle;     /* energy of saddle point */
  char          *structure;   /* structure and */
  float         energy;       /* .. energy of local min */
  float         Z;            /* partition function for structures in basin */
  float         Zg;           /* partition function for structures in basin */
  unsigned long my_GradPool;  /* for Gradient Basins */
  unsigned long my_pool;      /* # of structures attracted by this lmin */
  unsigned long fathers_pool; /* # of structures attracted by father-lmin */
  hash_entry    *left;
  hash_entry    *right;
  int           *POV;   /* only for the partial order */
  char          global; /* mark global Pareto points */
} loc_min;

typedef enum {
  Barriers_no_rates,
  Barriers_binary_rates,
  Barriers_text_rates,
  Barriers_both_rates = Barriers_binary_rates | Barriers_text_rates,
} barriers_rates_type;

typedef struct {
  int                 print_saddles;
  int                 bsize;
  int                 ssize;
  unsigned long       max_print;
  double              minh;
  double              kT;
  int                 want_quiet;
  int                 want_verbose;
  int                 want_connected;
  char                *GRAPH, *MOVESET;
  FILE                *INFILE;
  char                *seq;
  int                 poset;
  int                 label;
  barriers_rates_type rates;
  double              noLP_rate;
  int                 microrates;
  char                *text_rates_file;
  char                *binary_rates_file;
} barrier_options;

typedef struct {
  hash_entry  *hp;
  char        key[128];
  int       num;
} path_entry;

typedef struct {
  char          *structure;   /* unpacked structuture */
  unsigned long n;            /* index in energy sorted list */
  unsigned long min;          /* minimum */
  unsigned long truemin;      /* truemin */
  unsigned long gradmin;      /* gradient minimum */
  unsigned long truegradmin;  /* true gradient minimum */
  float         energy;       /* energy of structure */
  float         min_energy;
  float         truemin_energy;
  float         gradmin_energy;
  float         truegradmin_energy;
} map_struc;

#endif
