#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include "rcs_llist.h"


rcs_llist *rcs_llist_create(void) {
  rcs_llist *l;

  l = (rcs_llist *) malloc(sizeof(rcs_llist));
  assert(l);

  l->v_j = llist_create();
  l->r = llist_create();
  
  return l;
}

void rcs_llist_destroy(rcs_llist **l) {
  llist_node *n_ptr, *n_next_ptr;
  
  assert(*l);
  assert((*l)->r);

  if ((*l)->v_j != NULL) {
    llist_destroy(&(*l)->v_j);
  }
  
  n_ptr = (*l)->r->head;

  while (n_ptr != NULL) {
    n_next_ptr = n_ptr->next;
    free(n_ptr);
    n_ptr = n_next_ptr;
  }
  free((*l)->r);
  
  free(*l);
  *l = NULL;
}

void rcs_llist_add_to_row(const rcs_llist *l, float v, int j) {
  rcs_llist_node *n_ptr;
  
  assert(l);

  n_ptr = (rcs_llist_node *) malloc(sizeof(rcs_llist_node));
  assert(n_ptr);
  
  n_ptr->v = v;
  n_ptr->j = j;
  
  llist_append(l->v_j, (void *) n_ptr);
}

void rcs_llist_add_new_row(const rcs_llist *l, float v, int j) {
  rcs_llist_node *n_ptr;

  assert(l);

  rcs_llist_add_to_row(l, v, j);

  n_ptr = (rcs_llist_node *) l->v_j->tail;

  llist_append(l->r, (void *) n_ptr);
}

void rcs_llist_printf(const rcs_llist *l) {
  llist_node *r_ptr, *n_next_ptr, *n_current_ptr;
  rcs_llist_node *v_j_node_ptr;
  
  assert(l);

  r_ptr = (llist_node *) l->r->head;

  while (r_ptr != NULL) {
    n_current_ptr = r_ptr->data;

    if (r_ptr->next) {
      n_next_ptr = r_ptr->next->data;
    }
    else {
      n_next_ptr = NULL;
    }
      
    while(n_current_ptr != n_next_ptr) {
      v_j_node_ptr = (rcs_llist_node *) n_current_ptr->data;
      printf("(%f, %d)\t", v_j_node_ptr->v, v_j_node_ptr->j);
      n_current_ptr = n_current_ptr->next;
    }

    printf("\n");
    
    r_ptr = r_ptr->next;
  }
}

void rcs_llist_export_bin(const rcs_llist *l, FILE *fid_v, FILE *fid_j, FILE *fid_r, int *n_elem_exported) {
  llist_node *r_ptr, *n_next_ptr, *n_current_ptr;
  rcs_llist_node *v_j_node_ptr;
  int count, n_r, index, i;
    
  assert(l);
  
  r_ptr = (llist_node *) l->r->head;

  index = 0;
  
  while (r_ptr != NULL) {
    i = index + *n_elem_exported;
    count = fwrite(&i, sizeof(int), 1, fid_r);
    assert(count == 1);
    
    n_current_ptr = r_ptr->data;

    if (r_ptr->next) {
      n_next_ptr = r_ptr->next->data;
    }
    else {
      n_next_ptr = NULL;
    }

    n_r = 0;
    
    while(n_current_ptr != n_next_ptr) {
      v_j_node_ptr = (rcs_llist_node *) n_current_ptr->data;

      count = fwrite(&v_j_node_ptr->v, sizeof(float), 1, fid_v);
      assert(count == 1);

      count = fwrite(&v_j_node_ptr->j, sizeof(int), 1, fid_j);
      assert(count == 1);
      
      n_current_ptr = n_current_ptr->next;
      n_r++;
    }

    index += n_r;
    
    r_ptr = r_ptr->next;
  }

  i = index + *n_elem_exported;
  count = fwrite(&i, sizeof(int), 1, fid_r);
  assert(count == 1);

  *n_elem_exported += index;
}
