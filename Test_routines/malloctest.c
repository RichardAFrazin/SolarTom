/* gcc malloctest.c -o malloctest */

/*#include <malloc.h>*/
#include <stdio.h>
#include <stdlib.h>

int
main(int argc, char **argv)
{
  int i, b[2][3];
  char *p;

  i = 1;
  p = malloc(1048576);
  fprintf(stderr,"%d ",i);
  while (p != NULL){
    i++;
    p = malloc(1048576);
    fprintf(stderr,"%d ",i);
  }
  fprintf(stderr,"\n");
  return(0);
}

