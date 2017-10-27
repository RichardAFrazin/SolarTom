#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "llist.h"


void printf_string(void *data) {
  char *s;

  s = (char *) data;

  printf("%s\n", s);
}

int main(int argc, char **argv) {
  llist *l;
  char *s;

  
  l = llist_create();

  s = (char *) malloc(sizeof(char) * 16);
  assert(s);
  strcpy(s, "string #1");
  llist_append(l, (void *) s);

  s = (char *) malloc(sizeof(char) * 16);
  assert(s);
  strcpy(s, "string #2");
  llist_append(l, (void *) s);

  s = (char *) malloc(sizeof(char) * 16);
  assert(s);
  strcpy(s, "string #3");
  llist_append(l, (void *) s);
  
  llist_printf(l, &printf_string);
  
  llist_destroy(&l);
  
  return 0;
}
