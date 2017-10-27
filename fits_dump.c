/*  This is Butala's old eit_image.c.
 *    It read a .fts file and dumps the contents into unformatted binary.
 *
 * Mark D. Butala
 * Oct. 29, 2004
 *
 * Print CRPIX1, CRPIX2 and then dump image data to a file as an
 * unformatted string.
 *
 *  to compile:
 *gcc fits_dump.c -lm  -I/Users/frazin/tomography/wcstools-3.7.0/libwcs -L/Users/frazin/tomography/wcstools-3.7.0/libwcs -lwcs   -o fits_dump
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <fitsfile.h>

#define SIZE 1024

typedef float IM_TYPE;


void usage(char *command_name) {
  printf("usage: %s <fits file> <image output file> [size] \n", command_name);
}

int main(int argc, char **argv) {
  char *fits_file_name, *output_file_name;
  IM_TYPE *im_ptr;
  FILE *output_file;
  char *fits_head, *fits_image;
  int fits_lhead, fits_nbhead;
  int i, j;
  int size = SIZE;
  double center_x, center_y;
  
  char optstring[] = "fh";
  char opt;

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

  if ((argc - optind) < 2) {
    usage(argv[0]);
    return 1;
  }

  switch (argc - optind) {
  case 3:
    size = strtol(argv[optind+2], (char **)NULL, 10);
  case 2:
    fits_file_name = argv[optind];
    output_file_name = argv[optind+1];
    break;
  default:
    usage(argv[0]);
    return 1;
  }

  assert(fits_head = fitsrhead(fits_file_name, &fits_lhead, &fits_nbhead));
  assert(fits_image = fitsrimage(fits_file_name, fits_nbhead, fits_head));

  assert(hgetr8(fits_head, "CRPIX1", &center_x)); 
  assert(hgetr8(fits_head, "CRPIX2", &center_y));
	
	
  fprintf(stdout,"CRPIX1 = %g, CRPIX2 = %g\n", center_x, center_y);
  
  assert(output_file = fopen(output_file_name, "w"));

/*
  assert(fwrite(&center_x, sizeof(double), 1, output_file) == 1);
  assert(fwrite(&center_y, sizeof(double), 1, output_file) == 1);
  assert(fwrite(&exp_time, sizeof(double), 1, output_file) == 1);
  */
	
  im_ptr = (IM_TYPE *) fits_image;

  for (j = 0; j < size; j++ ) {
    for (i = size - 1; i >=  0; i--) {
      assert(fwrite(im_ptr + i * size + j, sizeof(IM_TYPE), 1, output_file) == 1);
    }
  }

  fclose(output_file);

  return 0;
}
