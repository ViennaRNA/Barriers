/* stapel.h */

#ifndef BARRIERS_STAPEL_H
#define BARRIERS_STAPEL_H

void
push(char *form);


char *
pop(void);


int
get_top(void);


void
ini_stapel(int size);


void
free_stapel(void);


void
reset_stapel(void);


#endif
