#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "headers.h"


void usage(char *argv0) {
  printf("%s <-h> [<infile> <outfile> <matrix suffix> <lamnda1> <lambda2> <lambda3>]\n", argv0);
  printf("Either all or none of the optional parameters must be specified\n");
  printf("-v: Use verbose mode on callsolves\n");
}

int main(int argc, char **argv) {
  char command[MAXPATH];
  
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

  if (argc - optind == 0) {
    sprintf(command, "./callsolve_cg -v");
    printf("Executing: %s\n", command);
    fflush(stdout);
    assert(system(command) == 0);
    
    sprintf(command, "./callsolve_fess -v -0");
    printf("Executing: %s\n", command);
    fflush(stdout);
    assert(system(command) == 0);
  }
  else if (argc - optind == 6) {
    /* infile => outfile */
    sprintf(command, "./callsolve_cg -v %s %s %s %s %s %s",
            argv[optind], argv[optind+1], argv[optind+2],
            argv[optind+3], argv[optind+4], argv[optind+5]);
    printf("Executing: %s\n", command);
    fflush(stdout);
    assert(system(command) == 0);

    /* outfile => outfile */
    sprintf(command, "./callsolve_fess -v %s %s %s %s %s %s",
            argv[optind+1], argv[optind+1], argv[optind+2],
            argv[optind+3], argv[optind+4], argv[optind+5]);
    printf("Executing: %s\n", command);
    fflush(stdout);
    assert(system(command) == 0);
  }
  else {
    usage(argv[0]);
    return 1;
  }

  return 0;
}
