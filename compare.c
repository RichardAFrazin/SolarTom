/*  compare.c  rewritten by Richard Frazin summer 2008
 *
 *   This goes through a list of list of .fts files (set by 
 * CONFSTRING in buildA_params.h) and uses the spacecraft
 * coordinates derived from each to create synthetic intensities
 * as determined from a 3D reconstruction (in an input file).
 *
 *  Changes to include WISPRI/O by Alberto Vásquez, Fall 2017
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <fitsfile.h>
#include "headers.h"

#define IMAGE_DIR TOMROOT"bindata/Compare/"
#define COMPARISON_CALCULATION /*this is used by buildrow.c */
#define QEPS 1.e-4 /*this is used by buildrow.c */


char *fitsrhead(char *, int *, int *);
char *fitsrimage(char *, int, char *);


static void usage(char *arg0) {
  printf("usage: %s <x_infile> [k]\n", arg0);
  printf(" x_infile is the filename of reconstruction.\n");
  printf("[k]: Compute the k'th comparison (first is k=0)\n");
}


int main(int argc, char **argv){
  char confstring[] = CONFSTRING;
  FILE *fid_conf, *fid;
  char idstring[MAXPATH], filename[MAXPATH], BpBcode[] = "xx", filename_x[MAXPATH];
  const double rmax = RMAX;
  static float pBval[IMSIZE][IMSIZE], pBcalc[IMSIZE][IMSIZE];
  static float rho[IMSIZE][IMSIZE], eta[IMSIZE][IMSIZE];
  double sun_ob1[3], sun_ob2[3], spol1[3], sob[3], spol2[3], r3tmp[3];
  double dsun, pang, deltagrid, rho1, eta1, carlong, mjd, roll_offset;
  double dsun_obs, ddat, J2k_OBS[3], obslat, sun_ob3[3];// Extra variables added by Albert, mainly for testing purposes, but also sun_ob3 serves to determine sign of t3, the "time" of the spacecraft.
  Rot R12, R23, Rtmp, *Rx, *Ry, *Rz;
  int nfiles, i, ii, jj, kk, ll, mmm, hasdata, totalB;
  const int nc3 = NBINS, imsize = IMSIZE;
  float *Arow_long, *x;
  double center_x, center_y;
  const double pixsize = PIXSIZE;
  int k = -1, start, stop;
  int opt;
  char optstring[] = "h";

  x = (float *) malloc(NBINS * sizeof(float));
  Arow_long = (float *) malloc(NBINS * sizeof(float));
  if (x == NULL || Arow_long == NULL){
    fprintf(stderr,"compare.c: malloc error, x and/or Arow_long.\n");
    fflush(stderr);
  }

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

  if ((argc > 3) || (argc < 2)){
    usage(argv[0]);
    return(1);
  }
  
  strcpy(filename_x,argv[1]);
  if (argc == 3)
    k = atoi(argv[2]);

  /* load x vector */
  strcpy(filename, BINDIR);
  strcat(filename, filename_x);
  fid = fopen(filename, "rb");
  if (fid == NULL) {
    fprintf(stderr, "%s\n", filename);
    fprintf(stderr, "Input file(s) not found\n");
    exit(2);
  }
  if (fread(x, sizeof(float), nc3, fid) != nc3) {
    fprintf(stderr, "compare: problem reading file: %s\n",filename_x);
    exit(1);
  }
  fclose(fid);

  /* open the list file */
  if ((fid_conf = fopen(confstring, "r")) == NULL) {
    fprintf(stderr, "Main: Input file '%s' not found\n", confstring);
    exit(1);
  }
  fgets(idstring, MAXPATH, fid_conf);
  sscanf(idstring, "%d", &nfiles);
  fprintf(stderr, "There are %d files.\n", nfiles);

  if (k == -1) {
    start = 0;
    stop = nfiles;
  }
  else {
    if (k < 0) {
      fprintf(stderr, "k >= 0!\n");
      return 1;
    }
    if (k >= nfiles) {
      fprintf(stderr, "k < nfiles (%d)!\n", nfiles);
      return 1;
    }
    start = k;
    stop = k + 1;
    printf("start: %d\tstop: %d\n", start, stop);
  }

  for (ii = 0; ii < nfiles; ii++) {
    fscanf(fid_conf, "%s\n", idstring);    
    if (ii >= start && ii < stop) {
      fprintf(stderr, "\ncurrent file: %s, %d of %d files\n", idstring,ii + 1, nfiles);
      strcpy(filename,DATADIR);
      strcat(filename,idstring);
   {   /* start of block */
    char *image;	/* FITS image */
    char *header;	/* FITS header */
    int lhead;		/* Maximum number of bytes in FITS header */
    int nbhead;		/* Actual number of bytes in FITS header */
    PB_IMTYPE *pbvector;
    float dist;
    FILE *fid_test;
    float x_image[IMSIZE], y_image[IMSIZE];
    
    /* test for existence of the data file */
    fid_test = NULL;
    fid_test = fopen(filename,"r");
    if (fid_test == NULL){
      fprintf(stderr,"compare.c: can't open file: %s\n",filename);
      fprintf(stderr,"Is there white space at the end of the filename?\n"); 
      fflush(stderr);
      exit(11);
    }
    fclose(fid_test);

    if ((header = fitsrhead(filename, &lhead, &nbhead)) != NULL) {
      if ((image = fitsrimage(filename, nbhead, header)) == NULL) {
        fprintf(stderr, "Cannot read FITS image %s\n", filename);
        free(header);
        return(1);
      }
    } else {
      fprintf(stderr, "Cannot read FITS file %s\n", filename);
      return(1);
    }

    pbvector = (PB_IMTYPE *) image;

   /* get the center pixel */

    assert(hgetr8(header, "CRPIX1", &center_x) || hgetr8(header, "XSUN", &center_x));
    assert(hgetr8(header, "CRPIX2", &center_y) || hgetr8(header, "YSUN", &center_y));
    center_x -= 1.; /* CRPIX(1,2), (XY)SUN use the start-at-1 convention */ 
    center_y -= 1.; 

    /* Get the roll angle offset b/c North may not be at the top 
     *   of the image */

#if (defined C2BUILD || defined C3BUILD)  /* INITANG1 is in the Marseilles files */
    assert(hgetr8(header,"CROTA1",&roll_offset) || hgetr8(header, "INITANG1", &roll_offset));
#elif (defined WISPRIBUILD || defined WISPROBUILD)
    assert(hgetr8(header,"CROTA2",&roll_offset));
    assert(hgetr8(header,"DSUN_OBS",&dsun_obs));
    assert(hgetr8(header,"J2kX_OBS",&ddat));   J2k_OBS[0]=ddat;
    assert(hgetr8(header,"J2kY_OBS",&ddat));   J2k_OBS[1]=ddat;
    assert(hgetr8(header,"J2kZ_OBS",&ddat));   J2k_OBS[2]=ddat;
    assert(hgetr8(header,"CRLT_OBS",&obslat));
#elif defined EITBUILD
    assert(hgetr8(header,"SC_ROLL",&roll_offset));
#elif (defined EUVIBUILD || defined CORBUILD || defined AIABUILD)
    assert(hgetr8(header,"CROTA2" ,&roll_offset));
#endif

#ifdef MARSEILLES /* we used to think CROTA1 = INITANG1 - 0.5  (.5 deg offset), but not anymore */ 
    roll_offset -= 0. ;  /* 0.5; */
#endif

    totalB = 0; /*this avoids unused variable warning*/
    
#if (defined WISPRIBUILD || defined WISPROBUILD)
    totalB = 1;
#endif

#ifdef C2BUILD  /*for C2, is it pB or total B ?*/
  strncpy(BpBcode, idstring + 3, 2); 
  if ( strcmp(BpBcode,"PB") == 0 )
    totalB = 0;
  else if ( strcmp(BpBcode,"BK") == 0 )
    totalB = 1;
  else {
    totalB = -1;
    fprintf(stderr,"Bad BpBcode: %s, idstring: %s\n",BpBcode, idstring); 
   exit(1);
  }
#endif

    /* Get the sun --> observer vector in J2000 GCI coords (km)
     * and the Carrington longitude (deg)
     */
    get_orbit(idstring, sun_ob1, &carlong, &mjd);
    carlong = carlong * M_PI / 180.0;
    dist = sun_ob1[0] * sun_ob1[0] + sun_ob1[1] * sun_ob1[1] + sun_ob1[2] * sun_ob1[2];
    dist = sqrt(dist);
 
    for (i = 0; i < imsize; i++)
    {
      x_image[i] = pixsize * ((float) i - center_x);
      y_image[i] = pixsize * ((float) i - center_y);
    }

      /* Albert's test printouts */
  fprintf(stderr,"Pixsize:            %3.10g ARCSEC\n",pixsize);
  fprintf(stderr,"Center X and Y:     %3.10g, %3.10g PX\n",center_x,center_y);
  fprintf(stderr,"Min and Max X:      %3.10g, %3.10g ARCSEC\n",x_image[0],x_image[imsize-1]);
  fprintf(stderr,"Min and Max angles: %3.10g, %3.10g DEG\n\n",x_image[0]/3600.,x_image[imsize-1]/3600.);
    
    /* when the observed image has solar north rotated clockwise from the top,
     *   CROTA (for EUVI) is positive.   I assume that the same convention
     *   holds for SC_ROLL and CROTA1 */
   
         for ( i = 0;  i < imsize;  i++) {
  //     for ( i = 0;  i < 1;  i++) {
#ifdef TESTBAND
         for (jj = imsize/2-BAND_WIDTH_PX/2; jj < imsize/2+BAND_WIDTH_PX/2; jj++) {
  //     for ( jj = 64;  jj < 65;  jj++) {
#else
      for (jj = 0; jj < imsize; jj++) {
#endif
        rho[i][jj] = (float) sqrt((double) (x_image[i]*x_image[i] + y_image[jj]*y_image[jj]));
        eta[i][jj] = (float) atan2((double) (-x_image[i]), (double) y_image[jj]) 
	           + (float) (roll_offset*0.017453292519943);
	//	  fprintf(stderr,"rho[%d][%d] = %3.10g deg\n",i,jj,rho[i][jj]/3600.);
	//	  fprintf(stderr,"eta[%d][%d] = %3.10g deg\n"   ,i,jj,eta[i][jj]*57.295778);
	//	  fprintf(stderr,"Impact = %3.10g RSUN\n",sin(ARCSECRAD * rho[i][jj]) * dist/RSUN);
      }
    }
   
     // exit(-1);
     
    for ( i = 0;  i < imsize;  i++) {
#ifdef TESTBAND
      for (jj = imsize/2-BAND_WIDTH_PX/2; jj < imsize/2+BAND_WIDTH_PX/2; jj++) {
#else
      for (jj = 0; jj < imsize; jj++) {
#endif
	pBval[i][jj] = (float) *(pbvector + imsize * jj + i) ;
    }
    }

       for (i = 0; i < imsize; i++) {
//     for ( i = 0;  i < 2;  i++) {
#ifdef TESTBAND
      for (jj = imsize/2-BAND_WIDTH_PX/2; jj < imsize/2+BAND_WIDTH_PX/2; jj++) {
#else
      for (jj = 0; jj < imsize; jj++) {
#endif
	/* Keep only data within certain radius range  */
#if (defined C2BUILD || defined C3BUILD || defined CORBUILD || defined WISPRIBUILD || defined WISPROBUILD)
	//OLD CODE BY RICH----------------------------------------------------------
	//        if (( tan(ARCSECRAD * rho[i][jj]) * dist > INSTR_RMAX * RSUN ) ||
	//            ( tan(ARCSECRAD * rho[i][jj]) * dist < INSTR_RMIN * RSUN )) {
        //REPLACED BY NEW CODE BY ALBERT:
        if (( sin(ARCSECRAD * rho[i][jj]) * dist > INSTR_RMAX * RSUN ) ||
            ( sin(ARCSECRAD * rho[i][jj]) * dist < INSTR_RMIN * RSUN )) {
       //--------------------------------------------------------------------------
          pBval[i][jj] = -999.0;
#elif (defined EITBUILD || defined EUVIBUILD || defined AIABUILD)
     // Similar corrections here:
     //	if ( (tan(ARCSECRAD * rho[i][jj]) * dist  > INSTR_RMAX * RSUN ) ||
        if ( (sin(ARCSECRAD * rho[i][jj]) * dist  > INSTR_RMAX * RSUN ) ||
#ifdef RING_REJECT
     //	     ((tan(ARCSECRAD * rho[i][jj]) * dist  > INNER_REJECT_RAD * RSUN) &&
     //	      (tan(ARCSECRAD * rho[i][jj]) * dist  < OUTER_REJECT_RAD * RSUN))  ){
       	     ((sin(ARCSECRAD * rho[i][jj]) * dist  > INNER_REJECT_RAD * RSUN) &&
       	      (sin(ARCSECRAD * rho[i][jj]) * dist  < OUTER_REJECT_RAD * RSUN))  ){
#else 
	        0 ){
#endif
          pBval[i][jj] = -999.0;
#endif
        } else {
#if (defined C2BUILD || defined C3BUILD || defined CORBUILD)
          /* the .79 factor is to convert from units of mean brightness
           *    to 1.e10*(center brightness)           */
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *=
#ifdef NRL
	       1.e10 * 0.79;
#elif (defined MARSEILLES)
  	       0.79; /* Marseilles scaling */ 
#endif
#if (defined WISPRIBUILD || defined WISPROBUILD)
          /* Add needed factor (if needed) once we decide the units of the synthetic images */
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *= 1.
#endif	       
#endif 

#ifdef DROP_NEG_PB
          if (pBval[i][jj] < 0)
            pBval[i][jj] = -999.0;
#endif
        }
     } /* jj loop over image pixels */
    } /* i loop over image pixels */  
      
    free(image);
    }   /* end of block */
     
    dsun = 0.0;
    for (i = 0; i < 3; i++) {
      sun_ob1[i] /= RSUN;
      dsun += sun_ob1[i] * sun_ob1[i];
    }
    dsun = sqrt(dsun);

    /* Albert's test printouts */
    fprintf(stderr,"Computed dsun: %g Rsun\n",dsun);
    fprintf(stderr,"Header's dsun: %g Rsun\n\n",dsun_obs/(RSUN*1.e3));

    /* solar pole vector */
    spol1[0] = cos(DELTApo) * cos(ALPHApo);
    spol1[1] = cos(DELTApo) * sin(ALPHApo);
    spol1[2] = sin(DELTApo);
 
    /* Calculate R12 matrix:  R12 = Rx(a3)Ry(a2)Rz(a1) */

    /* Zero y component of sun_ob */
    Rz = rotz(-atan2(sun_ob1[1], sun_ob1[0]));
    rotvmul(sob, Rz, sun_ob1);
    /* Zero z component of Rz * sunob */
    Ry = roty(-atan2(sob[2], sob[0]));
    /* Zero y component of spol */
    rotmul(&Rtmp, Ry, Rz);
    rotvmul(spol2, &Rtmp, spol1);
    Rx = rotx(atan2(spol2[1], spol2[2]));
    /*fprintf(stderr,"compare: remove the 0 from line 297!\n");fflush(stderr);*/
  
    rotmul(&R12, Rx, &Rtmp);
    rotvmul(sun_ob2, &R12, sun_ob1);
    rotvmul(spol2, &R12, spol1);

    free(Rx);
    free(Ry);
    free(Rz);

    /* Calculate R23 matrix */
    /* solar axes in frame 2 */
    pang = atan2(spol2[0], spol2[2]);
    Ry = roty(pang);
    Rz = rotz(carlong);		/* correct */
    rotmul(&R23, Rz, Ry);
    free(Rz);
    free(Ry);

    // Compute Sun_ob3:
    rotvmul(sun_ob3, &R23, sun_ob2);

#if (defined C2BUILD || defined C3BUILD || defined WISPRIBUILD || defined WISPROBUILD)
  if (totalB == 1)
    fprintf(stderr,"Total Brightness image.\n");
  else 
    fprintf(stderr,"Polarized Brightness image.\n");
#endif

fprintf(stderr, "polar angle: %g radians = %g deg\n",
          pang, pang * 180. / ((double) M_PI));
    fprintf(stderr, "Carrington longitude: %1.12g radians =  %3.9g deg\n",
          carlong, carlong * 180. / ((double) M_PI));
    fprintf(stderr, "spol1: [%g, %g, %g]\n", spol1[0], spol1[1], spol1[2]);
    fprintf(stderr, "spol2: [%g, %g, %g]\n", spol2[0], spol2[1], spol2[2]);
    rotvmul(r3tmp, &R23, spol2);
    fprintf(stderr, "spol3: [%g, %g, %g]\n", r3tmp[0], r3tmp[1], r3tmp[2]);
    fprintf(stderr, "COMPUTED sun_ob1:        [%3.10g, %3.10g, %3.10g]\n",
	  sun_ob1[0], sun_ob1[1], sun_ob1[2]);
    fprintf(stderr, "HEADER'S J2000 sun_obs:  [%3.10g, %3.10g, %3.10g]\n",
 	  J2k_OBS[0]/RSUN/1.e3, J2k_OBS[1]/RSUN/1.e3, J2k_OBS[2]/RSUN/1.e3);

    fprintf(stderr, "sun_ob2: [%g, %g, %g]\n", sun_ob2[0], sun_ob2[1],
          sun_ob2[2]);

    fprintf(stderr, "      Computed sun_ob3:  [%3.10g, %3.10g, %3.10g]\n\n",sun_ob3[0], sun_ob3[1], sun_ob3[2]);
    /* loop over image pixels */


#ifdef TESTBAND
      for (kk = imsize/2-BAND_WIDTH_PX/2; kk < imsize/2+BAND_WIDTH_PX/2; kk++) {
#else
      for (kk = 0; kk < imsize; kk++) {
#endif
      
      for (ll = 0; ll < imsize; ll++) {

	for (mmm = 0; mmm < nc3; mmm++)
          Arow_long[mmm] = 0.0;
	pBcalc[ll][kk] = -999.0;
         
	if (fabs(pBval[ll][kk] + 999.0) > QEPS) {
          rho1 = ((double) rho[ll][kk]) * ARCSECRAD;
          eta1 =  (double) eta[ll][kk];
          pBcalc[ll][kk] = 0.0;
#include "buildrow.c"
	} /* if pBval != -999 */
      }   /* ll loop over image pixels*/
    }    /* kk loop over image pixels*/


    /* write output files */
    strcpy(filename, IMAGE_DIR);
    strcat(filename, "comp_");
    strcat(filename, filename_x);
    strcat(filename, "_");
    strcat(filename, idstring);
    filename[strlen(filename) - 4] = '\0';
    strcat(filename, ".dat");
    fid = fopen(filename, "wb");
    fwrite(pBcalc, sizeof(float), IMSIZE*IMSIZE, fid);
    fclose(fid);
    fprintf(stderr, "output file = %s\n", filename);
      
    strcpy(filename, IMAGE_DIR);
    strcat(filename, "orig_");
    strcat(filename, idstring);
    filename[strlen(filename) - 4] = '\0';
    strcat(filename, ".dat");
    fid = fopen(filename, "wb");
    fwrite(pBval, sizeof(float), IMSIZE*IMSIZE, fid);
    fclose(fid);
    fprintf(stderr, "input file = %s\n", filename);
   
   } /* start/stop if statement */
  } /* ii  loop over files */

  fclose(fid_conf);
  return(0);
}
