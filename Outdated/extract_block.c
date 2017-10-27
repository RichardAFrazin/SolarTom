/* extract_block.c
 *
 * by: Mark D. Butala  11/04
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <getopt.h>

#include "headers.h"


void usage(char *arg0) {
  printf("usage: %s <-h> <conf_file> <matrix_name> <block>\n", arg0);
}

int main(int argc, char **argv) {
  char buffer[MAXPATH];
  char matrix_block_name[MAXPATH], matrix_block_c_name[MAXPATH];
  char *conf_file_name, *matrix_name;
  FILE *fid_conf, *info_fid;
  int block, n_blocks;
  struct sparse A;
  float *y = NULL, *delta = NULL;
  int *block_start, row_start, row_end;
  
  int opt;
  char optstring[] = "h";

  
  /* process command line arguments */
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'h':
      usage(argv[0]);
      return 0;
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  if ((argc - optind) == 3) {
    conf_file_name = argv[optind];
    matrix_name = argv[optind+1];
    block = strtol(argv[optind+2], (char **)NULL, 10);

    if (block < 0) {
      fprintf(stderr, "<block> must be nonnegative\n");
      return 1;
    }
  }
  else {
    usage(argv[0]);
    return 1;
  }

  /* Open conf file and find out the number of blocks  */
  strcat(strcpy(buffer, BINDIR), conf_file_name);
  fid_conf = fopen(buffer, "r");
  assert(fid_conf != NULL);
  fscanf(fid_conf, "%d", &n_blocks);
  fclose(fid_conf);

  assert(block <= n_blocks - 1);
  
  /* Load A matrix */
  printf("Loading A (%s) ... ", matrix_name);
  fflush(stdout);
  load_sparse(&A, matrix_name, BINDIR);
  printf("done\n");
  
  /* Load y and delta vectors */
  printf("Loading y (y%s) ... ", matrix_name);
  fflush(stdout);
  sprintf(buffer, "%sy%s", BINDIR, matrix_name);
  load_vf(&y, A.nf, buffer);
  printf("done\n");
  printf("Loadind delta (delta_%s) ... ", matrix_name);
  fflush(stdout);
  sprintf(buffer, "%sdelta_%s", BINDIR, matrix_name);
  load_vf(&delta, A.nf, buffer);
  printf("done\n");

  /* Determine the block start indices from the info file */
  block_start = (int *) malloc((n_blocks + 1)* sizeof(int));
  assert(block_start != NULL);
  
  strcat(strcat(strcpy(buffer, BINDIR),"info_"), matrix_name);
  info_fid = fopen(buffer, "r");
  assert(info_fid != NULL);
  assert(fread(block_start, sizeof(int), n_blocks+1, info_fid) == (n_blocks + 1));
  fclose(info_fid);


  row_start = block_start[block];
  row_end = block_start[block+1]-1;

  free(block_start);


  sprintf(matrix_block_name, "%s_block", matrix_name);
  sprintf(matrix_block_c_name, "%s_block_c", matrix_name);

  printf("Block (small) matrix name: %s\n", matrix_block_name);
  printf("Block complement (large) matrix name: %s\n", matrix_block_c_name);
  
  printf("Extracting row %d-%d ... ", row_start, row_end);
  fflush(stdout);
  row_extract(A, y, delta, matrix_block_name, matrix_block_c_name, row_start, row_end);
  printf("done\n");
  
  free_vf(&y);
  free_vf(&delta);
  free_sparse(&A);
  
  return 0;
}
