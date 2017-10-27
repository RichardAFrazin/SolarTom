#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>


#define COL_BLOCK_SIZE 32768

static void convert_rcs_2_ccs(const char *rcs_v,
			      const char *rcs_j,
                              const char *rcs_r,
                              const char *ccs_v,
                              const char *ccs_i,
                              const char *ccs_n,
                              int nc3, int n_rows) {
  FILE *rcs_v_fid, *rcs_j_fid, *rcs_r_fid;
  FILE *ccs_v_fid, *ccs_i_fid, *ccs_n_fid;
  int done, j1, j2, i1, i2, i, k, n_elem, row_start, elem_count;
  float v_k;
  int j_k;
  int *i_ptr;
  float *v_ptr;
  llist_node *node_ptr;
  llist *v_list[COL_BLOCK_SIZE];
  llist *i_list[COL_BLOCK_SIZE];
  
  assert(rcs_v && rcs_j && rcs_r);
  assert(ccs_v && ccs_i && ccs_n);

  assert(rcs_v_fid = fopen(rcs_v, "r"));
  assert(rcs_j_fid = fopen(rcs_j, "r"));
  assert(rcs_r_fid = fopen(rcs_r, "r"));

  assert(ccs_v_fid = fopen(ccs_v, "w"));
  assert(ccs_i_fid = fopen(ccs_i, "w"));
  assert(ccs_n_fid = fopen(ccs_n, "w"));

  elem_count = 0;
  assert(fwrite(&elem_count, sizeof(int), 1, ccs_n_fid) == 1);
		
  done = 0;
  j1 = 0;
  j2 = j1 + COL_BLOCK_SIZE - 1;
  while (!done) {
    for (k = 0; k < COL_BLOCK_SIZE; k++) {
      v_list[k] = llist_create();
      i_list[k] = llist_create();
    }

    assert(fread(&i2, sizeof(int), 1, rcs_r_fid) == 1);

    for (i = 0; i < n_rows; i++) {
      i1 = i2;
      assert(fread(&i2, sizeof(int), 1, rcs_r_fid) == 1);
      n_elem = i2 - i1;
      row_start = i1;

      for (k = 0; k < n_elem; k++) {
	assert(fread(&v_k, sizeof(float), 1, rcs_v_fid) == 1);
	assert(fread(&j_k, sizeof(int), 1, rcs_j_fid) == 1);

	/* Is current element in the column block? */
	if (j_k >= j1 && j_k <= j2) {
	  assert(v_ptr = malloc(sizeof(float)));
	  assert(i_ptr = malloc(sizeof(int)));
	  *v_ptr = v_k;
	  *i_ptr = i;
	  llist_append(v_list[j_k-j1], (void *) v_ptr);
	  llist_append(i_list[j_k-j1], (void *) i_ptr);
	}
      }
    }

    for (k = 0; k < j2-j1+1; k++) {
      node_ptr = v_list[k]->head;
      while(node_ptr) {
	assert(fwrite(node_ptr->data, sizeof(float), 1, ccs_v_fid) == 1);
	node_ptr = node_ptr->next;
      }

      node_ptr = i_list[k]->head;
      while(node_ptr) {
	assert(fwrite(node_ptr->data, sizeof(int), 1, ccs_i_fid) == 1);
	node_ptr = node_ptr->next;
      }

      assert(v_list[k]->length == i_list[k]->length);
      elem_count += v_list[k]->length;
      assert(fwrite(&elem_count, sizeof(int), 1, ccs_n_fid) == 1);
    }
    
    for (k = 0; k < COL_BLOCK_SIZE; k++) {
      llist_destroy(&(v_list[k]));
      llist_destroy(&(i_list[k]));
    }

    rewind(rcs_v_fid);
    rewind(rcs_j_fid);
    rewind(rcs_r_fid);

    j1 = j1 + COL_BLOCK_SIZE;
    j2 = j2 + COL_BLOCK_SIZE;

    if (j2 >= nc3) {
      j2 = nc3 - 1;
    }
    
    done = (j1 >= nc3);
  }
  
  fclose(rcs_v_fid);
  fclose(rcs_j_fid);
  fclose(rcs_r_fid);

  fclose(ccs_v_fid);
  fclose(ccs_i_fid);
  fclose(ccs_n_fid);
}
