/* Simple program using Mink's library to extract data from a FITS file.
 *
 * Mark D. Butala
 * Jul. 2004
 *
 */



#include "headers.h"

#ifdef C2BUILD
#define C2
#elif defined C3BUILD
#define C3
#elif defined MK4BUILD
#define MK4
#elif defined AIABUILD
#define AIA
#else
#error No build target defined that is currently supported.
#endif


#if defined MK4 || defined MK4_POL || defined C3
typedef short IM_TYPE;
#endif

#ifdef C2
typedef float IM_TYPE;
#endif

#define BUFFER_SIZE 100

#define DATE_OBS_STRING "DATE-OBS"
#define TIME_OBS_STRING "TIME-OBS"

#define BSCALE_STRING   "BSCALE"
#define BZERO_STRING    "BZERO"

#define X_AXIS_STRING   "NAXIS1"
#define Y_AXIS_STRING   "NAXIS2"

#define THETA_AXIS_STRING       "NAXIS1"
#define R_AXIS_STRING           "NAXIS2"

#define R_SUN_STRING    "RSUN"

#define CENTER_X_STRING "CRPIX1"
#define CENTER_Y_STRING "CRPIX2"

#define INCR_X_STRING   "CDELT1"
#define INCR_Y_STRING   "CDELT2"

#define DISPMIN_STRING "DISPMIN"

#define CRRADIUS_STRING "CRRADIUS"

typedef unsigned short bool;

void usage(char *command_name) {
  printf("usage: %s [-v] [-m] [-s] <fits file> <image output file>\n", command_name);
  printf("-v: verbose mode\n");
  printf("-m: mask mode (MK4 and MK4_POL only)\n");
  printf("-s: scale mode (used with mask mode)\n");
}

int main(int argc, char **argv) {
  char optstring[] = "vms";
  char *fits_head;
  char *fits_file_name, *output_file_name;
  int fits_lhead, fits_nbhead;
  char buffer[BUFFER_SIZE], buffer2[BUFFER_SIZE];
  FILE *output_file;
  IM_TYPE *im_ptr;
  int year, month, day, hour, min, sec;
#if defined MK4 || defined MK4_POL

  int i, j;
  double bscale, bzero;
  float *scaled_im_ptr;
#endif
#if defined C2 || defined C3 || defined MK4

  short size_x, size_y;
#endif
#ifdef MK4_POL
  short size_r, size_theta;
  double center_r, incr_r, dispmin;
#endif

#ifdef MK4
  double center_x, center_y, r_sun, incr_x, incr_y, dispmin;
#endif
  
  char opt;
  bool verbose_mode, mask_mode, scale_mode;

  char *fits_image;

  verbose_mode = 0;
  mask_mode = 0;
  scale_mode = 0;
  
  while ((opt = getopt(argc, argv, optstring)) != -1) {
    switch (opt) {
    case 'v':
      verbose_mode = 1;
      break;
    case 'm':
#if defined C2 || defined C3
      fprintf(stderr, "Mask mode only available for MK4 images!\n");
      return 1;
#endif
      mask_mode = 1;
      break;
    case 's':
      scale_mode = 1;
      break;
    default:
      usage(argv[0]);
      return 1;
    }
  }

  if ((argc - optind) != 2) {
    usage(argv[0]);
    return 1;
  }

  fits_file_name = argv[optind];
  output_file_name = argv[optind+1];


  if (scale_mode && !mask_mode) {
    fprintf(stderr, "Scale mode requires mask mode!\n");
    return 2;
  }
  
  /* PROCESS HEADER */
  assert(fits_head = fitsrhead(fits_file_name, &fits_lhead, &fits_nbhead));

  assert(hgets(fits_head, DATE_OBS_STRING, BUFFER_SIZE, buffer));

#if defined C2 || defined C3
  /* extract year, month, and day from the header */
  year = strtol(strtok(buffer, "/"), NULL, 10);
  month = strtol(strtok(NULL, "/"), NULL, 10);
  day = strtol(strtok(NULL, "/"), NULL, 10);

  /* make sure date string is properly formated */
  assert(strtok(NULL, "/") == NULL);
#endif

#if defined MK4 || defined MK4_POL
  /* extract year, month, and day from the header */
  year = strtol(strtok(buffer, "-"), NULL, 10);
  month = strtol(strtok(NULL, "-"), NULL, 10);
  day = strtol(strtok(NULL, "-"), NULL, 10);

  /* make sure date string is properly formated */
  assert(strtok(NULL, "-") == NULL);
#endif


  assert(hgets(fits_head, TIME_OBS_STRING, BUFFER_SIZE, buffer));

#if defined C2 || defined C3
  /* extract the hour, minute, and second from the header */
  hour = strtol(strtok(buffer, ":"), NULL, 10);
  min = strtol(strtok(NULL, ":"), NULL, 10);
  sec = strtol(strtok(NULL, "."), NULL, 10);

  /* ignore the fraction of a second */
#endif

#if defined MK4 || defined MK4_POL
  /* extract the hour, minute, and second from the header */
  hour = strtol(strtok(buffer, ":"), NULL, 10);
  min = strtol(strtok(NULL, ":"), NULL, 10);
  sec = strtol(strtok(NULL, ":"), NULL, 10);

  /* make sure date string is properly formated */
  assert(strtok(NULL, ":") == NULL);
#endif

#if defined MK4 || defined MK4_POL
  /* extract the calibration parameters */
  assert(hgetr8(fits_head, BSCALE_STRING, &bscale));
  assert(hgetr8(fits_head, BZERO_STRING, &bzero));
#endif

#if defined C2 || defined C3 || defined MK4
  /* extract x and y dimensions of the data */
  assert(hgeti2(fits_head, X_AXIS_STRING, &size_x));
  assert(hgeti2(fits_head, Y_AXIS_STRING, &size_y));
#endif

#ifdef MK4_POL
  /* extract r and theta dimensions of the polar data */
  assert(hgeti2(fits_head, THETA_AXIS_STRING, &size_theta));
  assert(hgeti2(fits_head, R_AXIS_STRING, &size_r));
#endif


#ifdef MK4
  if (mask_mode) {
    assert(hgetr8(fits_head, CENTER_X_STRING, &center_x));
    assert(hgetr8(fits_head, CENTER_Y_STRING, &center_y));
    assert(hgetr8(fits_head, R_SUN_STRING, &r_sun));
    assert(hgetr8(fits_head, INCR_X_STRING, &incr_x));
    assert(hgetr8(fits_head, INCR_Y_STRING, &incr_y));
    assert(hgetr8(fits_head, DISPMIN_STRING, &dispmin));

    assert(incr_x == incr_y);
  }
#endif

#ifdef MK4_POL
  if (mask_mode) {
    assert(hgetr8(fits_head, CENTER_Y_STRING, &center_r));
    assert(hgetr8(fits_head, CRRADIUS_STRING, &incr_r));
    assert(hgetr8(fits_head, DISPMIN_STRING, &dispmin));
  }
#endif
  
  if (verbose_mode) {
    printf("dmy:\t\t%d\t%d\t%d\n", day, month, year);
    printf("hms:\t\t%d\t%d\t%d\n", hour, min, sec);

#if defined MK4 || defined MK4_POL

    printf("bzero:\t\t%g\n", bzero);
    printf("bscale:\t\t%g\n", bscale);
#endif

#if defined C2 || defined C3 || defined MK4

    printf("size_x:\t\t%d\n", size_x);
    printf("size_y:\t\t%d\n", size_y);
#endif

#ifdef MK4_POL
    printf("size_theta:\t%d\n", size_theta);
    printf("size_r:\t\t%d\n", size_r);
#endif

#ifdef MK4
    if (mask_mode) {
      printf("center_x:\t%f\n", center_x);
      printf("center_y:\t%f\n", center_y);
      printf("r_sun:\t\t%f\n", r_sun);
      printf("incr_x:\t\t%f\n", incr_x);
      printf("incr_y:\t\t%f\n", incr_y);
      printf("dispmin:\t%g\n", dispmin);
    }
#endif

#ifdef MK4_POL
    if (mask_mode) {
      printf("center_r:\t%f\n", center_r);
      printf("incr_r:\t\t%f\n", incr_r);
      printf("dispmin:\t%g\n", dispmin);
    }
#endif
  }


  /* OUTPUT IMAGE */
  strcpy(buffer, output_file_name);
#if defined C2 || defined C3 || defined MK4

  sprintf(buffer2, "_%dx%d", size_x, size_y);
#elif defined MK4_POL

  sprintf(buffer2, "_%dx%d", size_theta, size_r);
#endif

  output_file_name = strcat(buffer, buffer2);

  assert(fits_image = fitsrimage(fits_file_name, fits_nbhead, fits_head));

  assert(output_file = fopen(output_file_name, "w"));

  im_ptr = (IM_TYPE *) fits_image;


#ifdef C2

  /* Note that, strictly speaking, a transpose is required (as is done
   * below for MK4) so that the output of fitstest will display
   * oriented correctly in matlab after a reshape. */
  
  assert(fwrite(im_ptr, sizeof(float), size_x * size_y, output_file) == size_x * size_y);

#elif defined C3
  assert(fwrite(im_ptr, sizeof(short), size_x * size_y, output_file) == size_x * size_y);
  
#elif defined MK4
  
  assert(scaled_im_ptr = (float *) malloc(sizeof(float) * size_x * size_y));

  for (i = 0; i < size_y; i++) {
    for (j = 0; j < size_x; j++) {
      if (mask_mode) {
        double dist;

        dist = sqrt((j+1 - center_x)*(j+1 - center_x) + (i+1 - center_y)*(i+1 - center_y));


        /* Assumption: incr_x == incr_y.  This was checked with an
           assert when the values where obtained from the header. */

        if ((dist < (MK4_RMIN * r_sun / incr_x)) ||
            (dist > (MK4_RMAX * r_sun / incr_x))) {
          if (scale_mode) {
            scaled_im_ptr[i*size_x + j] = -999;
          }
          else {
            scaled_im_ptr[i*size_x + j] = dispmin;
          }
        }
        else {
          if (scale_mode) {
            scaled_im_ptr[i*size_x + j] = (im_ptr[j*size_y + i] * bscale + bzero) * 1e10 * 0.79;
          }
          else {
            scaled_im_ptr[i*size_x + j] = im_ptr[j*size_y + i] * bscale + bzero;
          }
        }
        
      }
      else {
        scaled_im_ptr[i*size_x + j] = im_ptr[j*size_y + i] * bscale + bzero;
      }
    }
  }

  assert(fwrite(scaled_im_ptr, sizeof(float), size_x * size_y, output_file) == size_x * size_y);

  free(scaled_im_ptr);
#elif defined MK4_POL

  assert(scaled_im_ptr = (float *) malloc(sizeof(float) * size_r * size_theta));


  for (i = 0; i < size_theta; i++) {
    for (j = 0; j < size_r; j++) {
      if (mask_mode) {
        double dist;
        
        dist = ((j+1) - center_r) / incr_r;

        if (dist < MK4_RMIN || dist > MK4_RMAX) {
          if (scale_mode) {
            scaled_im_ptr[i*size_r + j] = -999;
          }
          else {
            scaled_im_ptr[i*size_r + j] = dispmin;
          }
        }
        else {
          if (scale_mode) {
            scaled_im_ptr[i*size_r + j] = (im_ptr[j*size_theta + i] * bscale + bzero) * 1e10 * 0.79;
          }
          else {
            scaled_im_ptr[i*size_r + j] = im_ptr[j*size_theta + i] * bscale + bzero;
          }
        }
      }
      else {
        scaled_im_ptr[i*size_r + j] = im_ptr[j*size_theta + i] * bscale + bzero;
      }
    }
  }
  
  assert(fwrite(scaled_im_ptr, sizeof(float), size_r * size_theta, output_file) == size_r * size_theta);

  free(scaled_im_ptr);
#endif

  fclose(output_file);

  return 0;
}
