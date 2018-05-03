#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "headers.h"

static void squeeze(char *fname_d, char *fname_dd, int ncol, int rmin, int rmax);
static void v_split(float *v, int size_v,
                    int rmin, int rmax,
                    char *prefix, char *fname_d, char *fname_dd);
static void squeeze_mod(char *fname_d, char *fname_dd, int ncol, int rmin, int rmax, int rstart, int rend);
static void v_split_mod(float *v, int size_v,
		    int rmin, int rmax, int rstart, int rend,
                    char *prefix, char *fname_d, char *fname_dd);

/* vmultSparse calculates
 *   c = (A * b)*scalefactor
 *   where A is a matrix and b is a
 *   column vector.
 */

void
/*
vmultSparse(struct sparse a, float *b, float *c, float scalefactor)
*/
vmultSparse(a, b, c, scalefactor)
struct sparse a;
float *b, *c;
float scalefactor;
{
  float val;
  int i, k, indx;

  for (i = 0; i < a.nf; i++)
    *(c + i) = 0.0;

  for (i = 0; i < a.ncol; i++) {
    for (k = a.n[i]; k < a.n[i + 1]; k++) {
      indx = *(a.iB + k);
      val = *(a.vB + k);
      *(c + indx) += (*(b + i)) * val;
    }
  }

  for (i = 0; i < a.nf; i++)
    *(c + i) *= scalefactor;
}

/* vmult_tranposeSparse ADDS a'*b to c
 *   c must be initialized before call
 *   scalefactor has an obvious meaning
 *   0 <= begindex <= nc3-1 is the row
 *    of c (or col of a) corresponding
 *    to the beginning of the chunk. set
 *    to 0 if there is only 1 chunk.
 */

void vmult_transposeSparse(a, b, c, scalefactor, begindex)
struct sparse a;
float *b, *c;
float scalefactor;
int begindex;
{
  float val;
  int i, k, indx;

  for (i = 0; i < a.ncol; i++) {
    for (k = a.n[i]; k < a.n[i + 1]; k++) {
      indx = *(a.iB + k);
      val = *(a.vB + k);
      *(c + begindex + i) += *(b + indx) * val * scalefactor;
    }
  }
}

/* matmultSparse  c = a * b where
 *   a and b are suitably dimensioned
 *   matrices. c is not held in memory
 *   and is written to disk.  fileid
 *   contains the filename suffix.
 *   the prefixes are n, i, v
 */

void matmultSparse(a, b, fileid)
struct sparse a;
struct sparse b;
char *fileid;
{
  FILE *fid_nc, *fid_ic, *fid_vc;
  char fs1[MAXPATH], fs2[MAXPATH], fs3[MAXPATH], fs4[MAXPATH];
  int i, j, k;
  int indb, inda, count, max_inda, min_inda;
  float valb;
  float *val_vector;
  int *ind_vector;

  if (b.nf != a.ncol) {
    fprintf(stderr, "MATMULTSPARSE: misdimensioned matrices!\n");
    
    fprintf(stderr, "a.nf = %d\ta.ncol = %d\n", a.nf, a.ncol);
    fprintf(stderr, "b.nf = %d\tb.ncol = %d\n", b.nf, b.ncol);
    exit(69);
  }

  strcpy(fs4, fileid);
  strcpy(fs1, BINDIR);
  strcpy(fs2, BINDIR);
  strcpy(fs3, BINDIR);
  strcat(fs1, "n");
  strcat(fs1, fs4);
  strcat(fs2, "i");
  strcat(fs2, fs4);
  strcat(fs3, "v");
  strcat(fs3, fs4);
  fid_nc = fopen(fs1, "wb");
  fid_ic = fopen(fs2, "wb");
  fid_vc = fopen(fs3, "wb");
  if (fid_nc == NULL || fid_ic == NULL || fid_vc == NULL) {
    fprintf(stderr, "MATMULTSPARSE: file open error!\n files: ");
    fprintf(stderr, "%s\n%s\n%s\n", fs1, fs2, fs3);
    exit(0);
  }

  val_vector = malloc(a.nf * sizeof(float));
  ind_vector = malloc(a.nf * sizeof(int));
  if (val_vector == NULL || ind_vector == NULL) {
    fprintf(stderr, "MATMULTSPARSE: malloc error!\n");
    exit(1);
  }

  for (k = 0; k < a.nf; k++) {
    *(val_vector + k) = 0.0;
    *(ind_vector + k) = 0;
  }

  count = 0;
  fwrite(&count, sizeof(int), 1, fid_nc);

  for (i = 0; i < b.ncol; i++) {

    max_inda = 0;
    min_inda = a.nf - 1;
    if (b.n[i] == b.n[i + 1])
      min_inda = 0;

    for (j = b.n[i]; j < b.n[i + 1]; j++) {

      indb = b.iB[j];
      valb = b.vB[j];

      for (k = a.n[indb]; k < a.n[indb + 1]; k++) {
        inda = a.iB[k];
        if (inda > max_inda)
          max_inda = inda;
        if (inda < min_inda)
          min_inda = inda;
        *(val_vector + inda) += a.vB[k] * valb;
        *(ind_vector + inda) = 1;
      }
    }

    for (k = min_inda; k <= max_inda; k++) {
      if (*(ind_vector + k) != 0) {
        count += 1;
        fwrite(&k, sizeof(int), 1, fid_ic);
        fwrite(val_vector + k, sizeof(float), 1, fid_vc);
        *(val_vector + k) = 0.0;
        *(ind_vector + k) = 0;
      }
    }
    fwrite(&count, sizeof(int), 1, fid_nc);
  }


  fclose(fid_nc);
  fclose(fid_ic);
  fclose(fid_vc);

}


/* this is version of row_extract (see below). 
 *   rmin is the beginning row number to be removed.
 *   rmax        ending
 *   rstart is the beginning row to be written 
 *   rend            ending
 * note: rmin <= rstart < rend <= rmax
 */
void row_extract_mod(struct sparse A,
                 float *y, float *delta,
		 char *fname_d, char *fname_dd,
		  int rmin, int rmax, int rstart, int rend)
{
  int i, j, k;
  int *diag_d, *diag_dd;
  struct sparse d, dd;
  
  assert(rmin <= rstart && rstart <= rend && rend <= rmax);
  assert(rmax < A.nf);


/* The d matrix contains the rows from rmin to rmax.
   The dd matrix contains all the other rows. */

  /* set up diag vectors */
  diag_d = (int *) malloc(A.nf * sizeof(int));
  diag_dd = (int *) malloc(A.nf * sizeof(int));

  assert(diag_d != NULL && diag_dd != NULL);

  for (k = 0; k < A.nf; k++) {
    if(k >= rmin && k <= rmax) {
      diag_dd[k] = 0;
      diag_d[k]  = 0;
      if (k >= rstart && k <= rend)
	diag_d[k] = 1;
    }
    else {
      diag_d[k] = 0;
      diag_dd[k] = 1;
    }
  }

  /* set up d and dd in sparse structure */

  d.ncol = A.nf;
  dd.ncol = A.nf;
  d.nf = A.nf;
  dd.nf = A.nf;

  d.n = (int *) malloc((d.ncol + 1) * sizeof(int));
  dd.n = (int *) malloc((dd.ncol + 1) * sizeof(int));

  assert(d.n != NULL && dd.n != NULL);

  i = 0;
  for (k = 0; k < d.ncol; k++)
    if (diag_d[k] == 1)
      i++;

  j = 0;
  for (k = 0; k < dd.ncol; k++)
    if (diag_dd[k] == 1)
      j++;

  d.iB = (int *) malloc(i * sizeof(int));
  dd.iB = (int *) malloc(j * sizeof(int));

  assert(d.iB != NULL && dd.iB != NULL);
  
  d.vB = malloc(i * sizeof(float));
  dd.vB = malloc(j * sizeof(float));

  assert(d.vB != NULL && dd.vB != NULL);

  dd.n[0] = 0;
  d.n[0] = 0;

  j = 0;
  for (k = 0; k < d.ncol; k++){
    if (diag_d[k] == 1){
      d.n[k + 1] = d.n[k] + 1;
      d.iB[j] = k;
      d.vB[j] = 1;
      j++;
    } else {
      d.n[k + 1] = d.n[k];
    }
  }

  j = 0;
  for (k = 0; k < dd.ncol; k++){
    if (diag_dd[k] == 1){
      dd.n[k + 1] = dd.n[k] + 1;
      dd.iB[j] = k;
      dd.vB[j] = 1;
      j++;
    } else {
      dd.n[k + 1] = dd.n[k];
    }
  }
  

  /* multiply the matrices and write to disk */

  matmultSparse(d , A, fname_d);
  matmultSparse(dd, A, fname_dd);

  /* squeeze (remove rows of zeros) from the resultant matricies */

  squeeze_mod(fname_d, fname_dd, A.ncol, rmin, rmax, rstart, rend);

  /* split y and delta vectors */
  v_split_mod(y    , A.nf, rmin, rmax, rstart, rend, "y"     , fname_d, fname_dd);
  v_split_mod(delta, A.nf, rmin, rmax, rstart, rend, "delta_", fname_d, fname_dd);
  
  free(diag_d);
  free(diag_dd);
  free(d.n);
  free(d.iB);
  free(d.vB);
  free(dd.n);
  free(dd.iB);
  free(dd.vB);
}

 /* This code has been copied more or less verbatim from
  row_extract.c.  This could be done more efficiently, e.g. the diag
  vectors are not really needed, but the code has been debugged and
  works.  The squeeze code could be more efficient as well by
  incorporating it into matmultSparse, but it was easier this way. */

void row_extract(struct sparse A,
                 float *y,
                 float *delta,
                 char *fname_d, char *fname_dd,
                 int rmin, int rmax) {
  
  int i, j, k;
  int *diag_d, *diag_dd;
  struct sparse d, dd;
    
  assert(rmax >= rmin && rmin >= 0);
  assert(rmax < A.nf);

/* The d matrix contains the rows from rmin to rmax.
   The dd matrix contains all the other rows. */

  /* set up diag vectors */
  diag_d = (int *) malloc(A.nf * sizeof(int));
  diag_dd = (int *) malloc(A.nf * sizeof(int));

  assert(diag_d != NULL && diag_dd != NULL);

  for (k = 0; k < A.nf; k++) {
    if(k >= rmin && k <= rmax) {
      diag_d[k] = 1;
      diag_dd[k] = 0;
    }
    else {
      diag_d[k] = 0;
      diag_dd[k] = 1;
    }
  }

  /* set up d and dd in sparse structure */

  d.ncol = A.nf;
  dd.ncol = A.nf;
  d.nf = A.nf;
  dd.nf = A.nf;

  d.n = (int *) malloc((d.ncol + 1) * sizeof(int));
  dd.n = (int *) malloc((dd.ncol + 1) * sizeof(int));

  assert(d.n != NULL && dd.n != NULL);

  i = 0;
  j = 0;
  for (k = 0; k < d.ncol; k++) {
    if (diag_d[k] == 1) {
      i++;
    } else {
      j++;
    }
  }

  d.iB = (int *) malloc(i * sizeof(int));
  dd.iB = (int *) malloc(j * sizeof(int));

  assert(d.iB != NULL && dd.iB != NULL);
  
  d.vB = malloc(i * sizeof(float));
  dd.vB = malloc(j * sizeof(float));

  assert(d.vB != NULL && dd.vB != NULL);

  dd.n[0] = 0;
  d.n[0] = 0;

  i = 0;
  j = 0;
  for (k = 0; k < d.ncol; k++) {
    if (diag_d[k] == 1) {
      dd.n[k + 1] = dd.n[k];
      d.n[k + 1] = d.n[k] + 1;
      d.iB[j] = k;
      d.vB[j] = 1;
      j++;
    } else {
      d.n[k + 1] = d.n[k];
      dd.n[k + 1] = dd.n[k] + 1;
      dd.iB[i] = k;
      dd.vB[i] = 1;
      i++;
    }
  }

  /* multiply the matrices and write to disk */

  matmultSparse(d, A, fname_d);
  matmultSparse(dd, A, fname_dd);

  /* squeeze (remove rows of zeros) from the resultant matricies */

  squeeze(fname_d, fname_dd, A.ncol, rmin, rmax);

  /* split y and delta vectors */
  v_split(y    , A.nf, rmin, rmax, "y",      fname_d, fname_dd);
  v_split(delta, A.nf, rmin, rmax, "delta_", fname_d, fname_dd);
  
  free(diag_d);
  free(diag_dd);
  free(d.n);
  free(d.iB);
  free(d.vB);
  free(dd.n);
  free(dd.iB);
  free(dd.vB);
}


void squeeze(char *fname_d, char *fname_dd, int ncol, int rmin, int rmax) {
  FILE *d_n_fid, *dd_n_fid;
  FILE *d_i_fid, *dd_i_fid;
  char buffer[MAXPATH];
  int *n;
  int *i_d, *i_dd;
  int size_d, size_dd;
  int i, val;


  /* load d n vector */
  strcat(strcat(strcpy(buffer, BINDIR), "n"),fname_d);
  assert((d_n_fid = fopen(buffer, "r")) != NULL);
  assert((n = (int *) malloc(sizeof(int) * (ncol+1))) != NULL);
  assert(fread(n, sizeof(int), ncol+1, d_n_fid) == ncol+1);

  size_d = n[ncol];
  fclose(d_n_fid);

  /* load dd n vector */
  strcat(strcat(strcpy(buffer, BINDIR), "n"),fname_dd);
  assert((dd_n_fid = fopen(buffer, "r")) != NULL);
  assert(fread(n, sizeof(int), ncol+1, dd_n_fid) == ncol+1);

  size_dd = n[ncol];
  free(n);
  fclose(dd_n_fid);
  
  /* load d i vector for reading and writing */
  strcat(strcat(strcpy(buffer, BINDIR), "i"),fname_d);
  assert((d_i_fid = fopen(buffer, "r+")) != NULL);
  assert((i_d = (int *) malloc(sizeof(int) * size_d)) != NULL);
  assert(fread(i_d, sizeof(int), size_d, d_i_fid) == size_d);

  /* load dd i vector for reading and writing */
  strcat(strcat(strcpy(buffer, BINDIR), "i"),fname_dd);
  assert((dd_i_fid = fopen(buffer, "r+")) != NULL);
  assert((i_dd = (int *) malloc(sizeof(int) * size_dd)) != NULL);
  assert(fread(i_dd, sizeof(int), size_dd, dd_i_fid) == size_dd);

  rewind(d_i_fid);
  rewind(dd_i_fid);

  /* squeeze d */
  for (i=0; i<size_d; i++) {
    /* sanity check for d matrix */
    assert(i_d[i] >= rmin && i_d[i] <= rmax);
    
    val = i_d[i] - rmin;
    assert(fwrite(&val, sizeof(int), 1, d_i_fid) == 1);
  }

  /* squeeze dd */
  for (i=0; i<size_dd; i++) {
    /* sanity check for dd matrix */
    assert(i_dd[i] < rmin || i_dd[i] > rmax);
    
    if (i_dd[i] > rmax) {
      val = i_dd[i] - (rmax - rmin) - 1;
      assert(fwrite(&val, sizeof(int), 1, dd_i_fid) == 1);
    }
    else {
      val = i_dd[i];
      assert(fwrite(&val, sizeof(int), 1, dd_i_fid) == 1);
    }
  }

  /* cleanup */
  fclose(d_i_fid);
  fclose(dd_i_fid);
  free(i_d);
  free(i_dd);
}

void squeeze_mod(char *fname_d, char *fname_dd, int ncol, int rmin, int rmax, int rstart, int rend) {
  FILE *d_n_fid, *dd_n_fid;
  FILE *d_i_fid, *dd_i_fid;
  char buffer[MAXPATH];
  int *n;
  int *i_d, *i_dd;
  int size_d, size_dd;
  int i, val;


  /* load d n vector */
  strcat(strcat(strcpy(buffer, BINDIR), "n"),fname_d);
  assert((d_n_fid = fopen(buffer, "r")) != NULL);
  assert((n = (int *) malloc(sizeof(int) * (ncol+1))) != NULL);
  assert(fread(n, sizeof(int), ncol+1, d_n_fid) == ncol+1);

  size_d = n[ncol];
  fclose(d_n_fid);

  /* load dd n vector */
  strcat(strcat(strcpy(buffer, BINDIR), "n"),fname_dd);
  assert((dd_n_fid = fopen(buffer, "r")) != NULL);
  assert(fread(n, sizeof(int), ncol+1, dd_n_fid) == ncol+1);

  size_dd = n[ncol];
  free(n);
  fclose(dd_n_fid);
  
  /* load d i vector for reading and writing */
  strcat(strcat(strcpy(buffer, BINDIR), "i"),fname_d);
  assert((d_i_fid = fopen(buffer, "r+")) != NULL);
  assert((i_d = (int *) malloc(sizeof(int) * size_d)) != NULL);
  assert(fread(i_d, sizeof(int), size_d, d_i_fid) == size_d);

  /* load dd i vector for reading and writing */
  strcat(strcat(strcpy(buffer, BINDIR), "i"),fname_dd);
  assert((dd_i_fid = fopen(buffer, "r+")) != NULL);
  assert((i_dd = (int *) malloc(sizeof(int) * size_dd)) != NULL);
  assert(fread(i_dd, sizeof(int), size_dd, dd_i_fid) == size_dd);

  rewind(d_i_fid);
  rewind(dd_i_fid);

  /* squeeze d */
  for (i=0; i<size_d; i++) {
    /* sanity check for d matrix */
    assert(i_d[i] >= rstart && i_d[i] <= rend);
    
    val = i_d[i] - rstart;
    assert(fwrite(&val, sizeof(int), 1, d_i_fid) == 1);
  }

  /* squeeze dd */
  for (i=0; i<size_dd; i++) {
    /* sanity check for dd matrix */
    assert(i_dd[i] < rmin || i_dd[i] > rmax);
    
    if (i_dd[i] > rmax) {
      val = i_dd[i] - (rmax - rmin) - 1;
      assert(fwrite(&val, sizeof(int), 1, dd_i_fid) == 1);
    }
    else {
      val = i_dd[i];
      assert(fwrite(&val, sizeof(int), 1, dd_i_fid) == 1);
    }
  }

  /* cleanup */
  fclose(d_i_fid);
  fclose(dd_i_fid);
  free(i_d);
  free(i_dd);
}

void v_split(float *v, int size_v,
             int rmin, int rmax,
             char *prefix, char *fname_d, char *fname_dd) {
  char buffer[MAXPATH];
  int i;
  FILE *d_fid, *dd_fid;

  strcat(strcat(strcpy(buffer, BINDIR), prefix), fname_d);
  assert((d_fid = fopen(buffer, "w")) != NULL);

  strcat(strcat(strcpy(buffer, BINDIR), prefix), fname_dd);
  assert((dd_fid = fopen(buffer, "w")) != NULL);

  for (i=0; i < size_v; i++) {
    if (i >= rmin && i <= rmax) {
      assert(fwrite(&v[i], sizeof(float), 1, d_fid) == 1);
    }
    else {
      assert(fwrite(&v[i], sizeof(float), 1, dd_fid) == 1);
    }
  }

  fclose(d_fid);
  fclose(dd_fid);
}

void v_split_mod(float *v, int size_v,
	     int rmin, int rmax, int rstart, int rend,
             char *prefix, char *fname_d, char *fname_dd) {
  char buffer[MAXPATH];
  int i;
  FILE *d_fid, *dd_fid;

  strcat(strcat(strcpy(buffer, BINDIR), prefix), fname_d);
  assert((d_fid = fopen(buffer, "w")) != NULL);

  strcat(strcat(strcpy(buffer, BINDIR), prefix), fname_dd);
  assert((dd_fid = fopen(buffer, "w")) != NULL);

  for (i = 0; i < size_v; i++)
    if (i >= rstart && i <= rend)
      assert(fwrite(&v[i],sizeof (float), 1,  d_fid) == 1);

  for (i = 0; i < size_v; i++)
    if (i < rmin || i > rmax)
      assert(fwrite(&v[i],sizeof (float), 1, dd_fid) == 1);

  fclose(d_fid);
  fclose(dd_fid);
}


void load_sparse(struct sparse *A, char *name, char *dir) {
  FILE *fid_v, *fid_i, *fid_n;
  int k, max_row;
  static char buffer[MAXPATH];
  
  /* Load n vector */
  sprintf(buffer, "%sn%s", dir, name);

  /*assert((fid_n = fopen(buffer, "r")) != NULL);*/
  fid_n = fopen(buffer, "r");
  if (fid_n == NULL){
    fprintf(stderr,"\n sparse.c: load_sparse: Can't open %s\n",buffer);
    fflush(stderr);
    assert(0);
  }


  A->n = (int *) malloc(sizeof(int) * (NBINS+1));
  assert(fread(A->n, sizeof(int), NBINS+1, fid_n) == NBINS+1);


  A->ncol = NBINS;
  
  /* Load v vector */
  sprintf(buffer, "%sv%s", dir, name);
  assert((fid_v = fopen(buffer, "r")) != NULL);

  A->vB = (float *) malloc(sizeof(float) * A->n[NBINS]);

  assert(fread(A->vB, sizeof(float), A->n[NBINS], fid_v) == A->n[NBINS]);

  /* Load i vector */
  sprintf(buffer, "%si%s", dir, name);
  assert((fid_i = fopen(buffer, "r")) != NULL);

  A->iB = (int *) malloc(sizeof(int) * A->n[NBINS]);

  assert(fread(A->iB, sizeof(int), A->n[NBINS], fid_i) == A->n[NBINS]);
  
  
  /* Determine nf: # of rows */
  for (k = 0, max_row = 0; k < A->n[NBINS]; k++) {
    if (A->iB[k] > max_row) {
      max_row = A->iB[k];
    }
  }

  max_row++;

  A->nf = max_row;

  fclose(fid_n);
  fclose(fid_v);
  fclose(fid_i);
}

void free_sparse(struct sparse *A) {
  assert(A->n != NULL && A->vB != NULL && A->iB != NULL);

  free(A->n);
  free(A->vB);
  free(A->iB);

  A->n = A->iB = NULL;
  A->vB = NULL;
  A->nf = A->ncol = 0;
}

void load_vf(float **v, int n, char *fname) {
  FILE *fid;

  *v = (float *) malloc(sizeof(float) * n);
  assert(*v != NULL);

  fid = fopen(fname, "r");
  assert(fid != NULL);

  assert(fread(*v, sizeof(float), n, fid) == n);
  fclose(fid);
}

void free_vf(float **v) {
  assert(*v != NULL);
  free(*v);
  *v = NULL;
}
