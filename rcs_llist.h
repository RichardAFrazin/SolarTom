#ifndef RCS_LLIST_H
#define RCS_LLIST_H

#include "llist.h"


typedef struct {
  float v;
  int j;
} rcs_llist_node;

typedef struct {
  llist *v_j;
  llist *r;
} rcs_llist;


rcs_llist *rcs_llist_create(void);
void rcs_llist_destroy(rcs_llist **l);
void rcs_llist_add_to_row(const rcs_llist *l, float v, int j);
void rcs_llist_add_new_row(const rcs_llist *l, float v, int j);
void rcs_llist_printf(const rcs_llist *l);
void rcs_llist_export_bin(const rcs_llist *l, FILE *fid_v, FILE *fid_j, FILE *fid_r, int *n_elem_exported);

#endif
