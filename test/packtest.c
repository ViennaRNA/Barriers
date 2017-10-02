/* gcc packtest.c -o packtest ../utils.o */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"

  
  int main (void) {
    int i;

    while (1){
      char *line = get_line(stdin);
      if(line == NULL)
	break; /*EOF*/
      char *token=strtok(line," \t");
      if(token[0]=='A' || token[0]=='U' || token[0]=='G' || token[0]=='C')
	continue; /*skip sequence line*/
      char *s = strdup(token);
      int l = strlen(s);
      printf("%s<\n",s);

      char *key = pack_structure(s);
      /* printf("%s %s\n",s[i],key); */   
      char *struc = unpack_structure(key);
      printf("%s\n", struc);
   
      if (strncmp(s,struc,l) != 0 ){
	fprintf(stderr, "unmatched structures:\n>%s<\n>%s<\n",s,struc);
	exit(EXIT_FAILURE);
      }
      free(struc);
      free(key);
      free(s);
      printf("----\n");

    }

  return 0;
}
