#include <stdio.h>

#include "rcs_llist.h"


int main(int argc, char **argv) {
  rcs_llist *l1, *l2;

  l1 = rcs_llist_create();

  rcs_llist_add_new_row(l1, 1, 1);
  rcs_llist_add_to_row(l1, 2, 2);
  rcs_llist_add_to_row(l1, 3, 3);
  rcs_llist_add_new_row(l1, 4, 4);
  rcs_llist_add_new_row(l1, 5, 5);
  rcs_llist_add_to_row(l1, 6, 6);

  printf("l1:\n");
  rcs_llist_printf(l1);


  l2 = rcs_llist_create();

  rcs_llist_add_new_row(l2, 7, 7);
  rcs_llist_add_to_row(l2, 8, 8);
  rcs_llist_add_new_row(l2, 9, 9);
  rcs_llist_add_to_row(l2, 10, 10);

  printf("\nl2:\n");
  rcs_llist_printf(l2);

  /*
  rcs_llist_append(l1, l2);
  printf("\nl1 after append:\n");
  rcs_llist_printf(l1);
  */
  
  rcs_llist_destroy(&l1);
  
  return 0;
}
