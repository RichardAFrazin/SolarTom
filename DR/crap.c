#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>
#include <fitsfile.h>
#include "headers.h"

int main(int argc, char **argv) {
  char confstring[] = CONFSTRING;
  char a_outfile[MAXPATH];
  FILE *fid_y, *fid_delta, *fid_w, *fid_j, *fid_m;
  FILE *fid_conf, *fid_log, *fid_date, *fid_block;
  int nfiles, i, k, er;
  char idstring[MAXPATH], filestring[MAXPATH];
  char filename_v[MAXPATH], filename_i[MAXPATH], filename_n[MAXPATH];
  char fn_w[MAXPATH], fn_m[MAXPATH], fn_j[MAXPATH];
  int nc3, nA[NBINS], nA_old[NBINS], n_elem_exported;
  float *yy, *dd;
  int len_yy, len_y, int_m1, i1, *ip1;
  float float_m1, f1, *fp1;
  rcs_llist *rcs;
  int opt;
  char optstring[] = "h";
  double mjdi, mjd_av; /* for DR computation */
  char filename[MAXPATH], *fitsdate, *header;
  FILE *fid_test;
  int lhead, nbhead; 
  char example_date[] = "2009-03-20T15:01:00.000";
  mjd_av = -1.0;

  strcpy(a_outfile, A_OUTFILE);
  er = 8;
 
  return (0);
}
