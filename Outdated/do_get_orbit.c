#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <getopt.h>

#include "headers.h"

void usage(char *arg0) {
  printf("usage: %s [year] [month] [day] [hour] [minute] <second>\n", arg0);
}

int main(int argc, char **argv) {
  int year, month, day, hour, min, sec;
  double carlong, sun_ob[3];
  
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

  if (argc - optind <= 4) {
    usage(argv[0]);
    return 1;
  } else if ((argc - optind) <= 6) {
    sscanf(argv[optind], "%d", &year);
    sscanf(argv[optind+1], "%d", &month);
    sscanf(argv[optind+2], "%d", &day);
    sscanf(argv[optind+3], "%d", &hour);
    sscanf(argv[optind+4], "%d", &min);

    if ((argc - optind) == 6) {
      sscanf(argv[optind+5], "%d", &sec);
    }
    else {
      sec = 0;
    }
  }
  else {
    usage(argv[0]);
    return 1;
  }

  get_orbit_wrapper(year, month, day, hour, min, sun_ob, &carlong);
  printf("%+g\t%g\t%g\t%g\n",  carlong, sun_ob[0], sun_ob[1], sun_ob[2]);

  return 0;
}
