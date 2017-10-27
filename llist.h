#ifndef LLIST_H
#define LLIST_H

typedef struct llist_node_struct {
  void *data;
  struct llist_node_struct *next;
} llist_node;

typedef struct {
  llist_node *head;
  llist_node *tail;
  int *length;
} llist;

llist *llist_create(void);
void llist_destroy(llist **l);
void llist_append(llist *l, const void *data);
void llist_printf(const llist *l, void (*printf_node)(void *data));

#endif
