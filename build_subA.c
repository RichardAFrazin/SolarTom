/*
 *    build_subA:
 *  build the submatrix associated with each 2D image.
 * 
 *  by Paul Janzen and Richard Frazin Summer/Fall 1999
 *
 *  Changes to include WISPRI/O by Alberto Vásquez, Fall 2017
 *  Changes to handle LAM LASCO-C2 new headers by Alberto Vásquez, Fall 2017
 *  Changes to handle KCor, by Alberto Vásquez, February 2018
 *  Changes to handle CoMP, by Alberto Vásquez, May 2018
 *
 */

/* WARNING!!!
 *
 * Any change made in this file must also be made in compare.c!
 *
 */

/*#include <malloc.h>*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fitsfile.h>
#include <assert.h>

#include "headers.h"

#define QEPS 1.e-4
#ifndef M_PI
#define M_PI 3.14159265358979
#endif
 
/* pB and CoMP units may be expected to be different for each instrument, see bellow.
 * All of the ray computations are done in units of Rsun (6.96e5 km) 
 */

void build_subA(char *idstring, rcs_llist *rcs,
                float **y, float **delta, int *count, 
	        FILE *fid_log, FILE *fid_date) {
  const double rmax = RMAX;
  char filename[MAXPATH], BpBcode[] = "xx";
  static float pBval[IMSIZE][IMSIZE];
  static float rho[IMSIZE][IMSIZE];
  static float eta[IMSIZE][IMSIZE];
  float pB1, n_los;
  double sun_ob1[3], sun_ob2[3], spol1[3], sob[3], spol2[3], r3tmp[3];
  double dsun, pang, rho1, eta1, carlong, mjd, covar_factor, roll_offset;
  double dsun_obs, ddat, J2k_OBS[3], obslat, sun_ob3[3],tilt;// Extra variables added by Albert, mainly for testing purposes, but also sun_ob3 serves to determine sign of t3, the "time" of the spacecraft.
  Rot R12, R23, Rtmp;
  Rot *Rx, *Ry, *Rz;
  int i, jj, kk, ll, k, l, modnum, mmm, hasdata, yn, totalB;
  static float Arow_long[NBINS];// no-static causes seg fault /* i loop over columns */ near the end of this code
  const int nc3 = NBINS;
  const int imsize = IMSIZE;
  const int binfac = BINFAC;
  double center_x, center_y;
  const double pixsize = PIXSIZE;

  /* initialization stuff */

  totalB = 0;

#if (defined WISPRIBUILD || defined WISPROBUILD)
    totalB = 1;
#endif

#if (defined METISVLBUILD) // I am assuming only pB images from Metis.
    totalB = 0;
#endif

  strcpy(filename, DATADIR);
  strcat(filename, idstring);

#ifdef C2BUILD  /*for C2, is it pB or total B ?*/
#ifdef NRL
  strncpy(BpBcode, idstring + 3, 2); 
#endif
#ifdef MARSEILLES
  strncpy(BpBcode, idstring + 25, 2); // for new files
//strncpy(BpBcode, idstring + 20, 2); // for old files (CR-2208)
#endif
fprintf(stderr,"BpBcode: %s, idstring: %s\n",BpBcode, idstring);
  if ( strcmp(BpBcode,"PB") == 0 || strcmp(BpBcode,"pB") == 0)
    totalB = 0;
  else if ( strcmp(BpBcode,"BK") == 0 )
    totalB = 1;
  else {
    totalB = -1;
    fprintf(stderr,"Bad BpBcode: %s, idstring: %s\n",BpBcode, idstring); 
   exit(1);
  }
#endif

  yn = (imsize / binfac) * (imsize / binfac);
  *y = (float *) malloc(yn * sizeof(float));
  *delta = (float *) malloc(yn * sizeof(float));
  if (*y == NULL || *delta == NULL) {
    fprintf(stderr, "Malloc error in initialization\n");
    exit(1);
  }
  for (i = 0; i < yn; i++) {
    (*y)[i] = 0.0;
    (*delta)[i] = 0.0;
  }

  *count = 0;

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
      fprintf(stderr,"build_subA.c: can't open file: %s\n",filename);
      fprintf(stderr,"Is there white space at the end of the filename?\n"); 
      fflush(stderr);
      exit(11);
    }
    fclose(fid_test);

    if ((header = fitsrhead(filename, &lhead, &nbhead)) != NULL) {
      if ((image = fitsrimage(filename, nbhead, header)) == NULL) {
        fprintf(stderr, "Cannot read FITS image %s\n", filename);
        free(header);
        return;
      }
    } else {
      fprintf(stderr, "Cannot read FITS file %s\n", filename);
      return;
    }

    pbvector = (PB_IMTYPE *) image;

   /* get the center pixel, XSUN and YSUN are in the Marseilles files */

    assert(hgetr8(header, "CRPIX1", &center_x) || hgetr8(header, "XSUN", &center_x));
    assert(hgetr8(header, "CRPIX2", &center_y) || hgetr8(header, "YSUN", &center_y));
    center_x -= 1.; /* CRPIX(1,2), (XY)SUN use the start-at-1 convention */ 
    center_y -= 1.; 

    /* Get the roll angle offset (in deg) b/c North may not be at the top 
     *   of the image */

#if (defined C2BUILD || defined C3BUILD) /* INITANG1 used to be in the Marseilles files, not in their newest version. */
    assert(hgetr8(header,"CROTA1",&roll_offset) || hgetr8(header, "INITANG1", &roll_offset) ||  hgetr8(header, "ROLLANGL", &roll_offset));
  //assert(hgetr8(header,"R_SOHO",&dsun_obs)); THIS WAS THE KEYWORD FOR NRL LASCO C2, WE SWITCH TO THE NEXT MARSEILLE KEYWORD:
  //dsun_obs = dsun_obs*(RSUN*1.e3);// now in [m]
    assert(hgetr8(header,"DSUN"     ,&dsun_obs));    // [m]
    assert(hgetr8(header,"CRLT_OBS" ,&obslat));      // [deg]
#elif (defined WISPRIBUILD || defined WISPROBUILD)
    assert(hgetr8(header,"CROTA2"  ,&roll_offset));
    assert(hgetr8(header,"DSUN_OBS",&dsun_obs));
    assert(hgetr8(header,"J2kX_OBS",&ddat));   J2k_OBS[0]=ddat;
    assert(hgetr8(header,"J2kY_OBS",&ddat));   J2k_OBS[1]=ddat;
    assert(hgetr8(header,"J2kZ_OBS",&ddat));   J2k_OBS[2]=ddat;
    assert(hgetr8(header,"CRLT_OBS",&obslat));
#elif defined EITBUILD
    assert(hgetr8(header,"SC_ROLL",&roll_offset));
#elif (defined EUVIBUILD || defined CORBUILD || defined AIABUILD)
    assert(hgetr8(header,"CROTA2" ,&roll_offset));
    /*fprintf(stdout,"CROTA2 = %g\n",roll_offset); fflush(stdout);*/
#elif defined KCORBUILD
    assert(hgetr8(header,"INST_ROT" ,&roll_offset));
    if (roll_offset != 0.) {
      fprintf(stderr,"KCOR roll_offset = %g\n",roll_offset);
      exit(0);
    assert(hgetr8(header,"DSUN"     ,&dsun_obs));    // [m]
    assert(hgetr8(header,"CRLT_OBS" ,&obslat));      // [deg]
#elif defined METISVLBUILD
    assert(hgetr8(header,"INST_ROT" ,&roll_offset));
    if (roll_offset != 0.) {
      fprintf(stderr,"METIS roll_offset = %g\n",roll_offset);
      exit(0);
    }
    assert(hgetr8(header,"DSUN"     ,&dsun_obs));    // [m]
    assert(hgetr8(header,"CRLT_OBS" ,&obslat));      // [deg]
#elif defined COMPBUILD
    assert(hgetr8(header,"CROTA1"   ,&roll_offset));
    if (roll_offset != 0.) {
      fprintf(stderr,"CoMP roll_offset = %g\n",roll_offset);
      exit(0);
    }
    assert(hgetr8(header,"DSUN"     ,&dsun_obs));    // [m]
    assert(hgetr8(header,"CRLT_OBS" ,&obslat));      // [deg]
    /*
    assert(hgetr8(header,"CRLN_OBS" ,&carlong));     // [deg] 
    fprintf(stderr,"Header roll_offset: %e\n",roll_offset);
    fprintf(stderr,"Header DSUN:        %e\n",dsun_obs   );
    fprintf(stderr,"Header CRLT_OBS:    %e\n",obslat     );
    fprintf(stderr,"Header CRLN_OBS:    %e\n",carlong    );
    */
#endif

/*
#ifdef MARSEILLES // we used to think CROTA1 = INITANG1 - 0.5  (.5 deg offset), but not anymore
    roll_offset -= 0.5;
#endif
*/
       
    // Get the sun --> observer vector in J2000 GCI (CS-1) in [km]
    // and the sub-spacecraft Carrington longitude in [deg].
    get_orbit(idstring, sun_ob1, &carlong, &mjd);
    assert(fwrite(&mjd, sizeof(double), 1, fid_date) == 1);
    carlong = carlong * M_PI / 180.0;
    dist    = sun_ob1[0] * sun_ob1[0] + sun_ob1[1] * sun_ob1[1] + sun_ob1[2] * sun_ob1[2];
    dist    = sqrt(dist);

    /*  
    fprintf(stderr,"Computed dist: %3.10g Rsun\n",dist/RSUN);
    fprintf(stderr,"Header DSUN:   %3.10g Rsun\n",dsun_obs/RSUN/1.e3);    
    exit(0);
    */

    for (i = 0; i < imsize; i++)
    {
      x_image[i] = pixsize * ((float) i - center_x);
      y_image[i] = pixsize * ((float) i - center_y);
    }

    /* when the observed image has solar north rotated clockwise from the top,
       CROTA (for EUVI) is positive.   I assume that the same convention
       holds for SC_ROLL and CROTA1 */
   
    for ( i = 0;  i < imsize;  i++) {
    for (jj = 0; jj < imsize; jj++) {
        rho[i][jj] = (float) sqrt((double) (x_image[i]*x_image[i] + y_image[jj]*y_image[jj]));
        eta[i][jj] = (float) atan2((double) (-x_image[i]), (double) y_image[jj]) 
	           + (float) (roll_offset*0.017453292519943);	
    }
    }
    
    for ( i = 0;  i < imsize;  i++) {
    for (jj = 0; jj < imsize; jj++) {
	pBval[i][jj] = (float) *(pbvector + imsize * jj + i) ;
    }
    }

for (i = 0; i < imsize; i++) {
  for (jj = 0; jj < imsize; jj++) {
	/* Keep only data within certain radius range  */
#if (defined C2BUILD || defined C3BUILD || defined CORBUILD || defined WISPRIBUILD || defined WISPROBUILD || defined KCORBUILD || defined COMPBUILD || defined METISVLBUILD)
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
          /* the .79 factor is to convert from units of disk-mean brightness to disk-center brightness
           * the 1.e10 factor changes units from [Bsun] to [1E-10*Bsun] when appropriate (center brightness) */
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *=
#ifdef NRL
	      1.e10 * 0.79;    /* NRL scaling */ 
#elif (defined MARSEILLES)
  	       0.79;           /* Marseilles scaling */ 
#endif
#elif (defined WISPRIBUILD || defined WISPROBUILD)
          /* Add needed factor (if needed) once we decide the units of the synthetic images */
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *= 1.;
#elif (defined KCORBUILD)
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *= 1.e+10; // KCOR images are expected in units of [Bsun]
#elif (defined METISVLBUILD)
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *= 0.79; // METIS-VL images are expected in units of 1.e-10*<Bsun>
#elif (defined COMPBUILD)
	  if ( abs(pBval[i][jj] + 999) > QEPS)  /* check for -999 values (missing blocks) */
	    pBval[i][jj] *= 1.0; // Keep units of the data
#endif
	  
#ifdef DROP_NEG_PB
          if (pBval[i][jj] < 0)
            pBval[i][jj] = -999.0;
#endif
        }  /* if/else */
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
  fprintf(stderr,"Computed dsun: %3.10g Rsun\n",dsun);
  fprintf(stderr,"Header's dsun: %3.10g Rsun\n\n",dsun_obs/(RSUN*1.e3));

  // Solar pole vector in SC-1
  spol1[0] = cos(DELTApo) * cos(ALPHApo);
  spol1[1] = cos(DELTApo) * sin(ALPHApo);
  spol1[2] = sin(DELTApo);

  // Calculate R12 matrix as:  R12 = Rx(c) * Ry(b) * Rz(a)

  // Rz(a) will zero the "y" component of "sun_ob1"
  Rz = rotz(-atan2(sun_ob1[1], sun_ob1[0])); // Albert: "-atan2" within the "rotz" operator is correct,
                                             // as this is a clockwise rotation.
  rotvmul(sob, Rz, sun_ob1);

  fprintf(stderr, "            Sun_ob1: [%g, %g, %g]\n", sun_ob1[0], sun_ob1[1], sun_ob1[2]);
  fprintf(stderr, "      Rz(a) Sun_ob1: [%g, %g, %g]\n", sob[0],     sob[1],     sob[2]);

  // Ry(b) will zero the "z" component of "sob" ( = Rz(a) * sun_ob1 )
  Ry = roty(-atan2(sob[2], sob[0])); // Albert: The roty operator is defined in rots.c
                                     // to be clockwise for angle > 0, the sign "-"  here
                                     // is then correct, as this is a counter-clockwise rotation.

  rotmul(&Rtmp, Ry, Rz);

  rotvmul(sob, &Rtmp, sun_ob1);
  fprintf(stderr, "Ry(b) Rz(a) Sun_ob1: [%g, %g, %g]\n", sob[0],     sob[1],     sob[2]);  

  // R12 will zero the "y" component of spol1.
  // Also "R12 Sun_ob1" below should equal "Ry Rz Sun_ob1" above

  rotvmul(r3tmp, &Rtmp, spol1);         // Albert: Note this is NOT spol2 yet, it is: "Ry(b) Rz(a) spol1".
                                        // so I changed it to r3tmp here, for clarity
  Rx = rotx(atan2(r3tmp[1], r3tmp[2])); // Albert: "+atan2" within the "rotx" operator is correct,
                                        // as this is a counter-clockwise rotation.

  rotmul(&R12, Rx, &Rtmp);
  rotvmul(sun_ob2, &R12, sun_ob1);
  rotvmul(spol2, &R12, spol1);

  fprintf(stderr, "        R12 Sun_ob1: [%g, %g, %g]\n", sun_ob2[0], sun_ob2[1], sun_ob2[2]);
  fprintf(stderr, "              spol2: [%g, %g, %g]\n", spol2[0], spol2[1], spol2[2]);
  
  free(Rx);
  free(Ry);
  free(Rz);
  //-----------------R12 computed------------------------------------------------------------------------

  // If dealing with KCOR or COMP or LASCO-C2 (MARSEILLE) or METIS-VL data: R12, spol2 and sun_ob2 above are crap.
  // As R12 was only needed to compute spol2 and sun_ob2, we just forget about it,
  // and simply re-compute spol2 and sun_ob2 using the sub-Earth latitude and the,
  // Observer-Sun distance, which are both known from the header:
#if (defined KCORBUILD || defined COMPBUILD || defined C2BUILD || defined METISVLBUILD)
    tilt     = obslat*M_PI/180.0;
    spol2[0] = sin(tilt); // Note that tilt>0 implies North-pole towards Earth.
    spol2[1] = 0.;
    spol2[2] = cos(tilt);
  sun_ob2[0] = dsun_obs/1.e3/RSUN; // sun_ob2 must be in Rsun units.
  sun_ob2[1] = 0.;
  sun_ob2[2] = 0.;
  fprintf(stderr,"tilt:  %e\n",tilt);
#endif
  
  // Calculate R23 matrix as:  R23 = Rz(CarLong) * Ry(Tilt)
  /* solar axes in frame 2 */
  pang = atan2(spol2[0], spol2[2]); // Albert: This is the correct expression for the tilt angle of spol_2.

  // Zero "x" component of spol2
  Ry = roty(pang);                  // Albert: The roty operator is defined in rots.c
                                    // to be clockwise for angle > 0, the sign "+" here 
                                    // is then correct, as this is a clockwise rotation.
 
  Rz = rotz(carlong); /* correct */ // Albert: "+carlong" within the "rotx" operator is correct,              
                                    // as this is a counter-clockwise rotation.
                                    // I believe Rich put the "/* correct */" comment because
                                    // F&J (2000) says Rz(-Carlong), which is not correct.
  rotmul(&R23, Rz, Ry);

  rotvmul(r3tmp, &R23, spol2);
  fprintf(stderr, "              spol3: [%g, %g, %g]\n", r3tmp[0], r3tmp[1], r3tmp[2]);

  // Compute Sun_ob3:
  rotvmul(sun_ob3, &R23, sun_ob2);
  
  free(Rz);
  free(Ry);
  //-----------------R23 computed------------------------------------------------------------------------


#if (defined C2BUILD || defined C3BUILD || defined WISPRIBUILD || defined WISPROBUILD)
  if (totalB == 1)
    fprintf(stderr,"Total Brightness image.\n");
  else 
    fprintf(stderr,"Polarized Brightness image.\n");
#endif

  fprintf(stderr, "Polar angle: %g radians = %g deg\n", pang, pang * 180. / ((double) M_PI));
  fprintf(stderr, "     Header's Observed Latitude = %g deg\n", obslat);
  
  fprintf(stderr, "Carrington longitude: %1.12g radians =  %3.9g deg\n", carlong, carlong * 180. / ((double) M_PI));
  /* fprintf(stderr, "spol1: [%g, %g, %g]\n", spol1[0], spol1[1], spol1[2]);*/
  /* fprintf(stderr, "spol2: [%g, %g, %g]\n", spol2[0], spol2[1], spol2[2]);*/
  rotvmul(r3tmp, &R23, spol2);
  /* fprintf(stderr, "spol3: [%g, %g, %g]\n", r3tmp[0], r3tmp[1], r3tmp[2]);*/
  fprintf(stderr, "COMPUTED sun_ob1:        [%3.10g, %3.10g, %3.10g]\n",
	  sun_ob1[0], sun_ob1[1], sun_ob1[2]);
  fprintf(stderr, "HEADER'S J2000 sun_obs:  [%3.10g, %3.10g, %3.10g]\n",
 	  J2k_OBS[0]/RSUN/1.e3, J2k_OBS[1]/RSUN/1.e3, J2k_OBS[2]/RSUN/1.e3);
  fprintf(stderr, "      Computed sun_ob3:  [%3.10g, %3.10g, %3.10g]\n\n",sun_ob3[0], sun_ob3[1], sun_ob3[2]);

  fprintf(stderr, "Sub-Spacecraft Latitude  computed as ATAN(sun_ob3[2]/Sqrt{sun_ob3[0]^2+sun_ob3[1]^2)}: %3.10g deg\n",
	  atan2(sun_ob3[2],(double)sqrt(sun_ob3[0]*sun_ob3[0]+sun_ob3[1]*sun_ob3[1])) *180./((double)M_PI) );
  fprintf(stderr, "Sub-Spacecraft Longitude computed as ATAN{sun_ob3[1]/sun_ob3[0]}:                       %3.10g deg\n",
	  atan2(sun_ob3[1],sun_ob3[0]) *180./((double)M_PI)  + 360.);
  
//fprintf(stderr, "sun_ob2: [%g, %g, %g]\n",  sun_ob2[0], sun_ob2[1], sun_ob2[2]);

  
  fprintf(fid_log,"cl= %g deg, polar_ang= %g deg, so1=[%g, %g, %g] Rs\n", 
	  carlong*180./((double) M_PI), pang*180./((double) M_PI),
          sun_ob1[0], sun_ob1[1], sun_ob1[2]);

  /*  build Arow_long and y, delta:
   *  This part calculates the "matrix elements" for a given
   *  image pixel and then adds all of those elements associated
   *  with a single bin of image pixels together to make a row of AA.
   *  First it tests to see if the pixel contains any data by comparing to -999.
   */

  if ((imsize % binfac) == 0) {
    modnum = 0;
  } else {
    modnum = binfac;
  }

  for (kk = 0; kk < imsize - modnum; kk += binfac) {
    for (ll = 0; ll < imsize - modnum; ll += binfac) {

      // Make Arow_long ALL zeros to start.
      for (mmm = 0; mmm < nc3; mmm++)
        Arow_long[mmm] = 0.0;
      
      hasdata = 0;
      n_los = 0.0;

      /*--------Loop LOSs within BIN---------------------------------*/     
      for (k = 0; k < binfac; k++) {
        for (l = 0; l < binfac; l++) {
          if (ll + l > (imsize - 1) || kk + k > (imsize - 1)) {
            fprintf(stderr,
                    "build_subA: image pixel out of bounds error!\n");
            fprintf(stderr, "kk = %d, ll = %d  ", kk, ll);
            fprintf(stderr, "k  = %d,  l = %d\n",  k,  l);
            exit(23);
          }

          pB1 = pBval[ll + l][kk + k];
          if (fabs(pB1 + 999.0) > QEPS) {

            hasdata = 1;  /* this will be set to 0 by buildrow.c 
                           * if the LOS misses the grid */
            rho1 = ((double) rho[ll + l][kk + k]) * ARCSECRAD;
            eta1 =  (double) eta[ll + l][kk + k];

#include "buildrow.c"

            if (hasdata == 1){
               covar_factor = 1.0;
               n_los += 1.0;  // Add up the number of LOS within BIN hitting the grid
               (*y)[*count] += pB1 / covar_factor; 
               (*delta)[*count] += DELTA;
	    } /* if hasdata */
          }
        }		/* l bin loop */
      }			/* k bin loop */
      /*------------End Loop LOSs within BIN-------------------------*/

      /*---Do this if at least one LOS within BIN has hit the grid---*/
      if (hasdata == 1) {
        int done_1;        
        (*y)[*count] /= n_los; 
        // Now y[BIN] contains the average of the values of all LOSs within the BIN.
        /* make sparse array */
        jj = 0;
        done_1 = 0;
	// Arow_long used below was computed in buildrow.c above. It contains as many
	// elements as VOXELS are in the GRID. It contains NON-zero values (columns)
	// that are crossed by the BIN LOSs.
	// Arow_long[i] should contain the A-matrix element of column "i" and row "count".
	// Note that buildrow computes this element n_los times for the BIN, adding the results,
	// with each call. It is normalized by n_los below to get an averaged value for the BIN.
        for (i = 0; i < nc3; i++) {
          if (Arow_long[i] != 0.0) {
            if (done_1) {
              rcs_llist_add_to_row(rcs, Arow_long[i] / n_los, i);  // normalized to n_los.
            }
            else {
              done_1 = 1;
              rcs_llist_add_new_row(rcs, Arow_long[i] / n_los, i); // normalized to n_los.
            }
          } // if arow_long[i] is NON-zero, i.e. if voxel "i" has been hit by BIN "count".
        }  /* i loop over columns */
        (*count)++; // Counts the number of BIN that has hit the grid 
      }	  /* if hasdata == 1 */
      /*----End actions in case at least one LOS within BIN has hit the grid---*/

    }  /* ll */
  }    /* kk */

  return;
}
