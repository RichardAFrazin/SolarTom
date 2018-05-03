#include <stdlib.h>
#include <assert.h>

#include "llist.h"


llist *llist_create(void) {
  llist *l;

  l = (llist *) malloc(sizeof(llist));
  assert(l);
  l->head = NULL;
  l->tail = NULL;
  l->length = 0;
  
  return l;
}

void llist_destroy(llist **l) {
  llist_node *n, *next;
  
  assert(*l);

  n = (*l)->head;

  while (n != NULL) {
    free(n->data);
    next = n->next;
    free(n);
    n = next;
  }

  free(*l);
  *l = NULL;
}

void llist_append(llist *l, const void *data) {
  llist_node *n;
  
  assert(l);
  assert(data);

  n = (llist_node *) malloc(sizeof(llist_node));
  assert(n);
  n->data = (void *) data;
  n->next = NULL;
  
  if (l->tail == NULL) {
    l->head = l->tail = n;
  }
  else {
    l->tail->next = n;
    l->tail = n;
  }

  l->length++;
}

void llist_printf(const llist *l, void (*printf_node)(void *data)) {
  llist_node *n;
  
  assert(l);
  assert(printf_node);

  n = l->head;

  while (n != NULL) {
    printf_node(n->data);
    n = n->next;
  }
}
