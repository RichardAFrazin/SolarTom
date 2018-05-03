#include <stdio.h>
#include <assert.h>

#include "headers.h"

int main(int argc, char **argv) {
  FILE *fid_info;
  char buffer[MAXPATH];
  int i;
  int n_blocks = 12;
  int row;
  unsigned char scratch;
  
  strcat(strcat(strcpy(buffer, BINDIR), "info_"), A_OUTFILE);
  fid_info = fopen(buffer, "r");
  assert(fid_info != NULL);

  for (i=0; i < n_blocks+1; i++) {
    assert(fread(&row, sizeof(int), 1, fid_info) == 1);
    printf("%d\n", row);
  }

  scratch = fgetc(fid_info);
  printf("%c\n", scratch);
  scratch = fgetc(fid_info);
  printf("%c\n", scratch);
  scratch = fgetc(fid_info);
  printf("%c\n", scratch);
  scratch = fgetc(fid_info);
  printf("%c\n", scratch);
  
  fclose(fid_info);
}
