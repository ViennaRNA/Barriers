/* Last changed Time-stamp: <2001-03-06 17:52:42 stadler> */
/* stapel.h */

#ifndef _stapel_h
#define _stapel_h

extern void push(char *form);
extern char *pop(void);
extern int get_top(void);
extern void ini_stapel(int size);
extern void free_stapel(void);
extern void reset_stapel(void);

#endif
