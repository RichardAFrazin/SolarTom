/* CARTESIAN: (x,y,z), (bin 0 in each coord. is betwen -RMAX and RMAX(1 + 2/NCELLS)
 * CYLINDRICAL: (r, phi, z) (radial bin 0 is between r = 0 and RMAX/NRAD,
 *                           phi    bin 0 is between 0 and  2pi/NPHI,
 *                           z      bin 0 is betwen -RMAX and RMAX(1 + 2/NCELLS))
 * HOLLOW_SPHERE: (r, theta, phi)
 *    (radial bin 0 is between RMIN and (RMAX-RMIN)/NRAD,
 *     theta  bin 0 is between -pi/2 and (-pi/2 + pi/NTHETA),
 *     phi    bin 0 is same as CYLINDRICAL case)
 */

/* geomtest.c uses RAYDIAGNOSE definition */

/* Edited by Alberto M. Vásquez, CLASP Fall-2017, to handle WISPR,
 * added comments for documentation and did cosmetic edits for readability.
 */

// Strategy to deal with the possibility of spacecraft being within computaional grid,
// by A.M.Vasquez. We compute the SIGNED times of the LOS entering the computational
// grid (t1), leaving it (t2), this itself has been much expanded to deal now in all possible
// geometrical situations. We also compute now the signed time corresponding to the spacecraft
// location (t3). All of these for now only implemented in hollow_sphere geometry.
// The three times are incorported in the array t. After all crossings have also been added,
// the array t is sorted in ascending order and then truncated to the following ranges:
// [t1,t2] if the spacecraft is outside the grid, [t3,t2] if it is inside.

/***********  BEGIN BUILDROW ***********************/
{
  static double nrpt[3], g1[3], unit[3], los1[3], los2[3];
  static double t1, t2, t3, arclength, xx, yy, zz, impact, r, vdhA, vdhB, vdhC, vdhD;
  static double deltagrid, grideps, rayeps;
  static double junk[6], t[NBINS], rtmp, ttmp, gam, sgam, cgam, ptmp;
  static int binbin[6], facedex[2], jij, tdex, index[3], ardex, ontarget;
  static double abstrmin, abstrmax, *dtpr; // New variables added by Albert
  static int index0; // New variables added by Albert
  static char case_string[]="0";

  rayeps = 1.e-6;
#if defined CYLINDRICAL || defined HOLLOW_SPHERE
  int wrap, binrmin;
  double rr, phiphi;
#endif
#ifdef CARTESIAN
  deltagrid = (2 * rmax) / (double) NCELLS;
  grideps = 1.e-6*deltagrid;
#elif defined CYLINDRICAL
  deltagrid = (2 * rmax) / (double) NZ;
  grideps = 1.e-6*MIN((rmax / (double) NRAD), deltagrid);
#elif defined HOLLOW_SPHERE
  deltagrid = (rmax - ((double) RMIN)) / (double) NRAD;
  grideps = 1.e-6*deltagrid;
#endif

#ifdef RAYDIAGNOSE
  fprintf(stderr,"ENTERING BUILDROW: rho1 = %1.12g, eta1 = %1.12g\n",rho1, eta1);
  fflush(stderr);
#endif

  /* unit is the LOS unit vector
   * nrpt is the nearest point vector
   *
   * unit points AWAY from the observer, TOWARDS the Sun
   *   (recall, in coordinate system 2, the x-axis points
   *   TOWARD the observer) */

  // Note in CASE-3-D Unit MAY point AWAY from the observer and AWAY the Sun,
  // I will ASSUME ight now this is already possible in the current code
  // and simply not contemplated so, and hence commented so, by Rich.
  // If so, pointing towards or away the Sun is simply characterized by
  // r3dot(unit,sun_ob3) > 0 or < 0, respectively.

  /* this is correct if eta1 is the usual solar PA (in radians) */
  Rx = rotx(eta1);

  //-----------------NRPT & UNIT-------------------------------------
  // These two calculations of NRPT and UNIT in CS-2 correspond
  // to Eqs. (9) and (10) in Frazin & Janzen (2002), respectively.
  // A.M.V. corrected the expressions for r3tmp[0] and r3tmp[2].
  // The error was an incorrect multiplicative factor 1/cos(rho1).
  // While being of up to order 1e-6 for EUV instruments, and up to 1e-4
  // for LASCO-C2, it would have been uo to order 1E-1 to 1E0 for WISPR.

  r3tmp[0] = dsun*sin(rho1)*sin(rho1);
  r3tmp[1] = 0.0;
  r3tmp[2] = dsun*sin(rho1)*cos(rho1);
  rotvmul(nrpt, Rx, r3tmp);

  g1[0] = -cos(rho1);
  g1[1] =  0.0;
  g1[2] =  sin(rho1);
  rotvmul(unit, Rx, g1);

  free(Rx);

  // Compute NRPT and UNIT in SC-3
  rotvmul(r3tmp, &R23, nrpt);
  r3eq(nrpt, r3tmp);
  rotvmul(r3tmp, &R23, unit);
  r3eq(unit, r3tmp);
  // From now on NRPT and UNIT are given in CS-3, in [Rsun] units.

  // LOS' impact parameter = Norm(NRPT)
  impact = sqrt(r3dot(nrpt, nrpt)); // [Rsun] units

  // Compute t3, the SIGNED "time" of the spacecraft location
  t3 = sqrt(dsun*dsun - impact*impact);
  if (r3dot(unit,sun_ob3) < 0)
    t3 *= -1.;
	  
  //------------------------------------------------------
  /* Calculate t1,t2, the "times" where ray enters and leaves computation region
   * los1 and los2 mark where the LOS enters and leaves the computation area
   *
   * junk[] are the (signed) distances from nrpt where the LOS crosses
   *    the max and min values of the computation box for each of the 3
   *    coordinates
   *
   * If the LOS misses the computation grid, set ontarget = 0 and
   *     exit buildrow.c .
   *
   * If any of the endpoints specified by the junk vector
   *    hits the edge of the grid, set ontarget = 1
   */

// Compute t1,t2 for CARTESIAN COORDINATES,
// as well as facedex[0,1]: faces of entry/exit.
#ifdef CARTESIAN
  for (jij = 0; jij < 6; jij += 2)
  {
    if ( fabs(unit[jij/2]) > rayeps )
    {
      // If unit[a] NE 0, with a=0,1,2 then set junk[a] and junk[a+1] so that
      // unit[a] * junk[a  ] + nrpt[a] = - rmax
      // unit[a] * junk[a+1] + nrpt[a] = + rmax
      junk[jij    ] = (-rmax - nrpt[jij / 2]) / unit[jij / 2] ;
      junk[jij + 1] = (+rmax - nrpt[jij / 2]) / unit[jij / 2] ;
    }
    else
    {
      // If unit[a] ~ 0, with a=0,1,2 then set
      // junk[a] = junk[a+1] = 1.e12, a huge number, so that
      // neither junk[a] nor junk[a+1] will contribute to determine facedex[] below.
      junk[jij    ] = 1.e12;
      junk[jij + 1] = 1.e12;
    }
  }
  ontarget = 0;
  facedex[0] = -1; /* facedex[0] contains junk[] of the 1st intersection */
  facedex[1] = -1; //        [1]                        2nd
  for (jij = 0; jij < 6; jij++)
  {
    // vector_g1 = vector_nrpt + junk[jij] * vector_unit.
    g1[0] = nrpt[0] + junk[jij]*unit[0];
    g1[1] = nrpt[1] + junk[jij]*unit[1];
    g1[2] = nrpt[2] + junk[jij]*unit[2];

    if ( (fabs(g1[0]) < rmax + grideps) &&
         (fabs(g1[1]) < rmax + grideps) &&
         (fabs(g1[2]) < rmax + grideps)   )
      {
            ontarget = 1;
	    if ( facedex[0] < 0){
	        facedex[0] = jij;
            } else {
	        facedex[1] = jij;
            }
      #ifdef  RAYDIAGNOSE
      fprintf(stderr,"intersection with computation cube face %d\n",jij);
      #endif
      }
  }
  if (ontarget == 0)
  goto salida;

// Compute t1,t2 for CYLINDRICAL COORDINATES:
// as well as facedex[0,1]: faces of entry/exit.
#elif defined CYLINDRICAL
  vdhA = (unit[0]*unit[0] + unit[1]*unit[1]);
  gam = 2.0 * (unit[0] * nrpt[0] + unit[1] * nrpt[1]);
  rtmp = nrpt[0] * nrpt[0] + nrpt[1] * nrpt[1] - rmax * rmax;
  if ( vdhA > rayeps) {
    junk[0] =
      (-gam - sqrt(gam * gam - 4. * vdhA * rtmp)) / (2.*vdhA);
    junk[1] =
      (-gam + sqrt(gam * gam - 4. * vdhA * rtmp)) / (2.*vdhA);
  } else {
    junk[0] = 1.e12;
    junk[1] = 1.e12;
  }
  junk[2] =  1.e12;	/* these will not contribute */
  junk[3] =  1.e12;
  for (jij = 4; jij < 6; jij += 2) {
    if ( fabs(unit[jij/2]) > rayeps ) {
      junk[jij] = (-rmax - nrpt[jij / 2]) / unit[jij / 2];
      junk[jij + 1] = (rmax - nrpt[jij / 2]) / unit[jij / 2] ;
    } else {
      junk[jij] = 1.e12;
      junk[jij + 1] = 1.e12;
    }
  }
  ontarget = 0;
  facedex[0] = -1; /* face[0] contains junk[] of the 1st intersection */
  facedex[1] = -1; /*     [2]                        2nd  */
  for (jij = 0; jij < 6; jij++){
    g1[0] = nrpt[0] + junk[jij]*unit[0];
    g1[1] = nrpt[1] + junk[jij]*unit[1];
    g1[2] = nrpt[2] + junk[jij]*unit[2];
#ifdef  RAYDIAGNOSE
        fprintf(stderr,"%d: z = %1.10g, rr = %1.10g\n",jij, g1[2],sqrt(g1[0]*g1[0] + g1[1]*g1[1]) );
        fflush(stderr);
#endif
    if ( (fabs(g1[2]) < rmax + grideps) &&
         (sqrt(g1[0]*g1[0] + g1[1]*g1[1]) < rmax + grideps) ){
            ontarget = 1;
	    if ( facedex[0] < 0){
	        facedex[0] = jij;
            } else {
	        facedex[1] = jij;
            }
#ifdef  RAYDIAGNOSE
            fprintf(stderr,"cylinder boundary %d: %1.10g, %1.10g, %1.10g\n",
		    jij,g1[0],g1[1],g1[2]);
            fflush(stderr);
#endif
    }
  }
  if ( ontarget == 0)
  goto salida;

// Compute t1,t2 for SPHERICAL COORDINATES:
#elif defined HOLLOW_SPHERE

  ontarget = 1;
  if (impact > RMAX ){
    /* Is the LOS outside of computation sphere?  If so, just
         treat it has having no data (see build_subA.c)*/
    ontarget = 0;
    goto salida;
  }

   // Un-signed crossing times at radii RMIN and RMAX.
   abstrmin = sqrt(((double) RMIN)*((double) RMIN) - impact*impact);
   abstrmax = sqrt(((double) RMAX)*((double) RMAX) - impact*impact);

   // Compute SIGNED t1 and t2 for all possible on-target gemoetrical situations:
   
   if (dsun > ((double) RMAX)){ // Cases 1.
     strcpy(case_string,"1");
     if (impact > 1.0){         // Cases 1A, 1B.
       junk[0] = -abstrmax;
       junk[1] =  abstrmax;
       }
     if (impact <= 1.0){        // Case 1C.
       junk[0] = -abstrmax;
       junk[1] = -abstrmin;
       }
   } // Cases 1.

   if (dsun >= ((double) RMIN) && dsun <= ((double) RMAX)){ // Cases 2.
     strcpy(case_string,"2");
     if (impact > 1.0){              // Cases 2A, 2B, 2C, 2D.
       junk[0] = -abstrmax;
       junk[1] =  abstrmax;
       }
     if (impact <= 1.0){             // Cases 2E, 2F.
       if (r3dot(unit,sun_ob3) < 0){ // Case 2E.
       junk[0] = -abstrmax;
       junk[1] = -abstrmin;
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 2F.
       junk[0] =  abstrmin;
       junk[1] =  abstrmax;
       }
       }
   } // Cases 2.

   if (dsun < ((double) RMIN)){      // Cases 3.
     strcpy(case_string,"3");
     if (impact > 1.0){              // Cases 3A, 3B.
       junk[0] =  abstrmin;
       junk[1] =  abstrmax;
       }
     if (impact < 1.0){
       if (r3dot(unit,sun_ob3) < 0){ // Case 3C; set ontarget=1 (LOS hits Sun w/o intersecting grid)
         ontarget = 0;
         goto salida;
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 3D.
       junk[0] =  abstrmin;
       junk[1] =  abstrmax;
       }
       }
   } // Cases 3

   // Assign the SIGNED values to t1 and t2:
   t1 = junk[0];
   t2 = junk[1];

#endif

#if defined CARTESIAN || defined CYLINDRICAL
   if (junk[facedex[0]] < junk[facedex[1]] ){
     t1 = junk[facedex[0]];
     t2 = junk[facedex[1]];
   } else {
     t1 = junk[facedex[1]];
     t2 = junk[facedex[0]];
   }
   /*  I'M NOT SURE THIS DOES ANY GOOD...
   t1 = t1 + grideps;
   if (t2 > 0){
     t2 = t2 - grideps;
   } else {
     t2 = t2 + grideps;
   }
   */
#endif

// LOS endpoints position vector coordinates in CS-3.
// This is valid for all geometries.
   for (jij = 0; jij < 3; jij++)
   {
    los1[jij] = nrpt[jij] + t1*unit[jij];
    los2[jij] = nrpt[jij] + t2*unit[jij];
   }

   // At this point: los1, los2 and sun_ob3, are all given in SC-3.

   // Compute distances Spacecraft-los1 (dlos1) and Spacecraft-los2 (dlos2)
   // not needed for now, leave it here commented out in case we need them in the future
   // dlos1 = sqrt(r3dot(r3sub(los1,sun_ob3),r3sub(los1,sun_ob3)));
   // dlos2 = sqrt(r3dot(r3sub(los2,sun_ob3),r3sub(los2,sun_ob3)));

  /* put the bin number of LOS endpoints into binbin array -
       see CARTESIAN example for ordering */

#ifdef CARTESIAN

  binbin[0] = floor((los1[0] + rmax) / deltagrid);
  binbin[1] = floor((los2[0] + rmax) / deltagrid);
  binbin[2] = floor((los1[1] + rmax) / deltagrid);
  binbin[3] = floor((los2[1] + rmax) / deltagrid);
  binbin[4] = floor((los1[2] + rmax) / deltagrid);
  binbin[5] = floor((los2[2] + rmax) / deltagrid);

  /* sometime small errors put the end point just
   *  outside the box */
  for (jij = 0; jij < 6; jij++){
    if (binbin[jij] == NCELLS)
      binbin[jij] = NCELLS - 1;
  }

#elif defined CYLINDRICAL

  rtmp = sqrt(los1[0]*los1[0] + los1[1]*los1[1]);
  binbin[0] = floor(rtmp * ((double) NRAD) / rmax);
  rtmp = sqrt(los2[0]*los2[0] + los2[1]*los2[1]);
  binbin[1] = floor(rtmp * ((double) NRAD) / rmax);

  /* sometime small errors put the end point just
   *  outside the computional cylinder */
  if (binbin[0] == NRAD)
    binbin[0] = NRAD - 1;
  if (binbin[1] == NRAD)
    binbin[1] = NRAD - 1;


  ptmp = atan2(los1[1], los1[0]);
  if (ptmp < 0.0)
    ptmp += 2.0 * M_PI;
  binbin[2] = floor(ptmp * ((double) NPHI) / 2.0 / M_PI);
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry phi = %3.6g deg, ",ptmp*180./M_PI);
#endif

  rtmp = atan2(los2[1], los2[0]);
  if (rtmp < 0.0)
    rtmp += 2.0 * M_PI;
  binbin[3] = floor(rtmp * ((double) NPHI) / 2.0 / M_PI);
  /* set the wrap parameter */
  wrap = 0;
  if ( fabs(rtmp - ptmp) > M_PI )
     wrap = 1;
#ifdef RAYDIAGNOSE
  fprintf(stderr,"exit phi = %3.6g deg ==> wrap = %d\n", rtmp*180./M_PI,wrap);
#endif

  binbin[4] = floor((los1[2] + rmax) / deltagrid);
  binbin[5] = floor((los2[2] + rmax) / deltagrid);

  /* sometime small errors put the end point just
   *  outside the computional cylinder */
  if (binbin[4] == NZ)
    binbin[4] = NZ - 1;
  if (binbin[5] == NZ)
    binbin[5] = NZ - 1;

#elif defined HOLLOW_SPHERE
   /* r, theta (polar angle), phi (azimuthal angle) coordinate order */

  // The next code computes binbin[0] and binbin[1] for all possible cases 1, 2 and 3,
  // defined above when computing t1 and t2.;
   if (dsun > ((double) RMAX)){ // Cases 1.
     if (impact > 1.0){         // Cases 1A, 1B.
       binbin[0] = NRAD-1;
       binbin[1] = NRAD-1;
     }
     if (impact <= 1.0){        // Case 1C.
       binbin[0] = NRAD-1;
       binbin[1] = 0;
       }
   } // Cases 1.

   if (dsun >= ((double) RMIN) && dsun <= ((double) RMAX)){ // Cases 2.
     if (impact > 1.0){              // Cases 2A, 2B, 2C, 2D.
       binbin[0] = NRAD-1;
       binbin[1] = NRAD-1;
       }
     if (impact <= 1.0){             // Cases 2E, 2F.
       if (r3dot(unit,sun_ob3) < 0){ // Case 2E.
       binbin[0] = NRAD-1;
       binbin[1] = 0;
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 2F.
       binbin[0] = 0;
       binbin[1] = NRAD-1;
       }
       }
   } // Cases 2.

   if (dsun < ((double) RMIN)){      // Cases 3.
     if (impact > 1.0){              // Cases 3A, 3B.
       binbin[0] = 0;
       binbin[1] = NRAD-1;
       }
     if (impact < 1.0){
       if (r3dot(unit,sun_ob3) < 0){ // Case 3C; set ontarget=1 (LOS hits Sun w/o intersecting grid)
         ontarget = 0;
         goto salida; // it went to salida already when computing t1,t2, just included the line for completeness of cases. 
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 3D.
       binbin[0] = 0;
       binbin[1] = NRAD-1;
       }
       }
   } // Cases 3

  // Compute Latitude [rad] of vector los1 in CS-3:
  rtmp = atan( los1[2] / sqrt( los1[0]*los1[0] + los1[1]*los1[1]) );
  // Compute Theta bin of vector los1, taking Theta=0 in South pole, and Theta=Pi in North pole.
  // Works well only wih uniform theta grid. !! 
  binbin[2] =  floor( (rtmp + M_PI/2.)*((double) NTHETA)/ M_PI );
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry theta = %g deg, ",rtmp*180./M_PI);
#endif
  // Same Latitude/Theta computations for vector los2:
  rtmp = atan( los2[2] / sqrt( los2[0]*los2[0] + los2[1]*los2[1]) );
  binbin[3] =  floor( (rtmp + M_PI/2.)*((double) NTHETA)/ M_PI );
#ifdef RAYDIAGNOSE
  fprintf(stderr,"exit theta = %g deg\n",rtmp*180/M_PI, wrap);
#endif
  // Compute Longitude [rad] of vector los1 in CS-3, making sure is in range [0,2pi].
  rtmp = atan2(los1[1], los1[0]);
  if (rtmp < 0.) rtmp += 2.*M_PI;
  // Compute Phi bin of vector los1
  binbin[4] = floor(rtmp * ((double) NPHI)/2./M_PI);
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry phi = %g deg, ",rtmp*180./M_PI);
#endif
  // Same Longitude/Phi computations for vector los2:
  ptmp = atan2(los2[1], los2[0]);
  if (ptmp < 0.) ptmp += 2.*M_PI;
  binbin[5] = floor(ptmp * ((double) NPHI)/2./M_PI);
  wrap = 0;
  if ( fabs(rtmp - ptmp) > M_PI )
    wrap = 1; // If unsigned Longitude difference between los1 and los2 is larger than Pi.
#ifdef RAYDIAGNOSE
  fprintf(stderr,"exit phi = %g deg, ==> wrap = %d\n",
      ptmp*180/M_PI, wrap);
#endif

#endif

  /* find the "times" of bin crossings */

  // First two elements are entry and exit "times"
  t[0] = t1;
  t[1] = t2;
  tdex = 2;

  // In hollow_sphere geometry, third element is the Spacecraft "time"
#ifdef HOLLOW_SPHERE
  t[2] = t3;
  tdex++;
#endif

#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry distance t1 = %g, exit distance t2 = %g,\n",t1,t2);
  fprintf(stderr,"nrpt = %1.10g %1.10g %1.10g\n",nrpt[0],nrpt[1],nrpt[2]);
  fprintf(stderr,"unit = %1.10g %1.10g %1.10g\n",unit[0],unit[1],unit[2]);
  fprintf(stderr,"point of entry %1.10g, %1.10g, %1.10g\n",nrpt[0]+ t1*unit[0], nrpt[1]+t1*unit[1], nrpt[2]+t1*unit[2]);
  fprintf(stderr,"point of exit  %1.10g, %1.10g, %1.10g\n",nrpt[0]+ t2*unit[0], nrpt[1]+t2*unit[1], nrpt[2]+t2*unit[2]);
  fprintf(stderr,"entry/exit bins (binbin) = %d %d %d %d %d %d\n",binbin[0],binbin[1],binbin[2],binbin[3],binbin[4],binbin[5]);
#endif

#ifdef CYLINDRICAL

  /* binrmin = bin of minimum radius */
  ttmp = -(unit[0]*nrpt[0] + unit[1]*nrpt[1]) /
          (unit[0]*unit[0] + unit[1]*unit[1]);
  rtmp = sqrt((nrpt[0] + ttmp*unit[0])*(nrpt[0] + ttmp*unit[0]) +
              (nrpt[1] + ttmp*unit[1])*(nrpt[1] + ttmp*unit[1]));
  binrmin = floor(rtmp*((double) NRAD) / rmax );

  rtmp = nrpt[2] + ttmp*unit[2];

  /*radial bin crossings */
  vdhA = (unit[0]*unit[0] + unit[1]*unit[1]);
  gam = 2.*(unit[0]*nrpt[0] + unit[1]*nrpt[1]);
  ptmp = nrpt[0]*nrpt[0] + nrpt[1]*nrpt[1];

#ifdef RAYDIAGNOSE
  fprintf(stderr,"RADIAL Bins: binrmin= %d, bin crossings: ",binrmin);
  fflush(stderr);
#endif

  for (jij = MAX(binbin[0], binbin[1]) - 1; jij >= binrmin; jij--) {
    rtmp = rmax * ((double) (jij + 1) / (double) NRAD); /* = outer radius of jij bin */
    rtmp = ptmp - rtmp * rtmp;
    ttmp =
      (-gam - sqrt(gam*gam - 4.*vdhA*rtmp)) / (2.*vdhA);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d,%g)",jij,ttmp);
      fflush(stderr);
#endif
    }
    ttmp =
      (-gam + sqrt(gam*gam - 4.*vdhA*rtmp)) / (2.*vdhA);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d,%g)",jij,ttmp);
      fflush(stderr);
#endif
    }
  }

  /* angular bin crossings */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nPHI Bins: bin crossings: ");
  fflush(stderr);
#endif

  if (wrap == 0) {
    for (jij = MIN(binbin[2], binbin[3]);
         jij < MAX(binbin[2], binbin[3]); jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0]*ptmp) / (unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 0: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij+1)*360./ (double) NPHI);
        exit(32);
      }
    }
  } else {
    for (jij = MAX(binbin[2], binbin[3]); jij < NPHI; jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0]*ptmp) / (unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij+1)*360./ (double) NPHI);
        exit(32);
      }
    }
    for (jij = 0; jij < MIN(binbin[2], binbin[3]); jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0]*ptmp) / (unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,((jij+1)*360./ (double) NPHI));
        exit(32);
      }
    }
  }

  /* z bin crossings */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nZ Bins: bin crossings: ");
  fflush(stderr);
#endif
  for (jij = MIN(binbin[4], binbin[5]);
       jij < MAX(binbin[4], binbin[5]); jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[2]) / unit[2] ;
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d,%g)",jij,ttmp);
      fflush(stderr);
#endif
    } else {
      fprintf(stderr,"ttmp out of bounds (Z bin crossing)!");
      fprintf(stderr,"ttmp = %g, ptmp = %g, jij = %d, z = %g \n",ttmp,ptmp, jij, (jij+1)*deltagrid - rmax);
      fflush(stderr);
      exit(32);
    }
  }

#elif defined HOLLOW_SPHERE

  /*radial bin crossings */
  if (impact <= ((double) RMIN)){
     binrmin = 0;  /* binrmin = bin of minimum radius */
  } else {
     binrmin = rad_bin_number(impact);
     //     fprintf(stderr,"impact = %g, bin = %d.\n",impact,binrmin);
  }

#ifdef RAYDIAGNOSE
  fprintf(stderr,"RADIAL Bins: binrmin= %d, bin crossings: ",binrmin);
  fflush(stderr);
#endif
  
  /* -2 because of the bin numbering and the entry point into the last bin is already marked by t1 */
  // or by t2 (case 2F and all cases 3).
  for (jij = NRAD - 2; jij >= binrmin; jij--) {
        dtpr = rad_bin_boundaries(jij);  
        rtmp = *dtpr; // outer boundary of cell jij.
	// fprintf(stderr,"index = %d, r = %g, dr = %g\n",jij,(*dtpr+*(dtpr+1))/2.,(*dtpr-*(dtpr+1)));
        ttmp = - sqrt(rtmp*rtmp - impact*impact) ; // take here the NEGATIVE root.
        t[tdex] = ttmp;
	tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
	// In this new code we ALWAYS keep also the positive root.
	// The whole thing works in all cases (1,2,3) now,
	// as the array t is filtered out for any t < t3 below.
	// OLD CODE:	if ( impact > 1. ) { // take the POSITIVE root also
	  ttmp *= -1.;
  	  t[tdex] = ttmp;
	  tdex++;
#ifdef RAYDIAGNOSE
          fprintf(stderr,"(%d,%g)",jij,ttmp);
          fflush(stderr);
#endif
	  //OLD CODE: }
  }

  //  exit(-1);
  /* polar angle bin crossings
   *
   * This formulation does not distinguish between positive
   *   and negative angles.  The times are the same and only
   *   the times matter in the end.  Just loop over negative
   *   angles.
   */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nTheta Bins: bin crossings: ");
  fflush(stderr);
#endif

  for (jij = 0; jij < NTHETA/2 + 1; jij++) {
	 gam = tan( (jij+1)*M_PI/((double) NTHETA) - M_PI/2. );
	 gam = gam*gam;
	 vdhA = gam*(unit[0]*unit[0] + unit[1]*unit[1]) - unit[2]*unit[2];
	 vdhB = 2.*( gam*(unit[0]*nrpt[0] + unit[1]*nrpt[1]) - unit[2]*nrpt[2]);
	 sgam = gam*(nrpt[0]*nrpt[0] + nrpt[1]*nrpt[1]) - nrpt[2]*nrpt[2];

//	 ttmp =  (- vdhB - sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA);
	 ttmp =  GridDivision(- vdhB - sqrt(vdhB*vdhB - 4.*vdhA*sgam),2.*vdhA);
		 
	 if ((ttmp > t1) && (ttmp < t2)) {
	    t[tdex] = ttmp;
	    tdex++;
#ifdef RAYDIAGNOSE
            fprintf(stderr,"(%d, %g deg, %g)"
	      ,jij,(jij+1)*180./((double) NTHETA) - 90.,ttmp);
            fflush(stderr);
#endif
	 }

//       ttmp =  (- vdhB + sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA);
	 ttmp =  GridDivision(- vdhB + sqrt(vdhB*vdhB - 4.*vdhA*sgam),2.*vdhA);
	 if ((ttmp > t1) && (ttmp < t2)) {
	    t[tdex] = ttmp;
	    tdex++;
#ifdef RAYDIAGNOSE
            fprintf(stderr,"(%d, %g deg, %g)"
	      ,jij,(jij+1)*180./((double) NTHETA) - 90.,ttmp);
            fflush(stderr);
#endif
	 }
  }

  /* azimuthal bin crossings */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nPHI Bins: bin crossings: ");
  fflush(stderr);
#endif

  if (wrap == 0) {
    for (jij = MIN(binbin[4], binbin[5]);
         jij < MAX(binbin[4], binbin[5]); jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
//    ttmp = (nrpt[1] - nrpt[0]*ptmp) / (unit[0]*ptmp - unit[1]);
      ttmp = GridDivision(nrpt[1] - nrpt[0]*ptmp , unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 0: out of bounds!!\n");
	fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp,jij,(jij+1)*360./ (double) NPHI);
	exit(32);
      }
    }
  } else {
    for (jij = MAX(binbin[4], binbin[5]); jij < NPHI; jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
//    ttmp = (nrpt[1] - nrpt[0]*ptmp)/(unit[0]*ptmp - unit[1]);
      ttmp = GridDivision(nrpt[1] - nrpt[0]*ptmp , unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij+1)*360./ (double) NPHI);
	exit(32);
      }
    } /* jij loop */
    for (jij = 0; jij < MIN(binbin[4], binbin[5]); jij++) {
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
//    ttmp = (nrpt[1] - nrpt[0]*ptmp) / (unit[0]*ptmp - unit[1]);
      ttmp = GridDivision(nrpt[1] - nrpt[0]*ptmp , unit[0]*ptmp - unit[1]);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij+1)*360./ (double) NPHI);
       	exit(32);
      }
    } /* jij loop */
  }

#elif defined (CARTESIAN)

  /* x bin crossings */
  for (jij = MIN(binbin[0], binbin[1]); jij < MAX(binbin[0], binbin[1]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[0])/unit[0];
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }
  /* y bin crossings */
  for (jij = MIN(binbin[2], binbin[3]); jij < MAX(binbin[2], binbin[3]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[1])/unit[1];
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }
  /* z bin crossings */
  for (jij = MIN(binbin[4], binbin[5]); jij < MAX(binbin[4], binbin[5]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[2])/unit[2];
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }

#endif

  // At this point, tdex = # of voxels that are threaded by the LOS, for all geometries.
  // +1 (spacecraft location) for hollow_sphere.
  
  // Sort all the bin crossing "times", in ascending order, works for all geometries:
  qsort((void *) t, tdex, sizeof(double), (void *) &doublecompare);

  // Find index0 such that t[index0]=t3
  index0 = 0;
  while (t[index0] < t3) index0++;
  //  fprintf(stderr, "t[%d] = %g.  t3 = %g.  t[%d] = %g.\n",index0,t[index0],t3,index0+1,t[index0+1]);
  
  // If case_string NE "2" then index0 = index0 + 1, to exclude spacecraft location.
  if (strcmp(case_string,"1") == 0 || strcmp(case_string,"3") == 0) index0++;

  //  fprintf(stderr, "Case # %s. Set index0 = %d.\n",case_string,index0);
  
  // Redefine array t as the elements of the original array with index >= index0, and adjust tdex.
  //  fprintf(stderr, "Old t[0] = %g. Old tdex = %d. Old t[tdex-1] = %g.\n"  ,t[0],tdex,t[tdex-1]);
  for (jij=index0 ; jij < tdex; jij++) t[jij-index0]  = t[jij];
  tdex   = tdex - index0;
  // fprintf(stderr, "New t[0] = %g. New tdex = %d. New t[tdex-1] = %g.\n\n",t[0],tdex,t[tdex-1]);      

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\ntimes: ");
  for (jij = 0;jij < tdex;jij++)
    fprintf(stderr,"%g ",t[jij]);

  fprintf(stderr,"\nVoxel Numbers: ");
  fflush(stderr);
#endif

    /* calculate the matrix elements */

  for (jij = 1; jij < tdex; jij++) {    // Loop through all voxels threaded by the LOS.
    // ttmp and arclength, valid for all geometries:
    ttmp = 0.5 * (t[jij] + t[jij - 1]); // "time" of the voxel middle point (VMP) of the LOS.
    arclength = t[jij] - t[jij - 1];    // voxel arclength (VAL) of the LOS.
    //    fprintf(stderr, "t[%d] = %g. t[%d] = %g. arclength = %g.\n",jij-1,t[jij-1],jij,t[jij],arclength);
    if (arclength < 0.0) {
      fprintf(stderr, "BUILDROW: arclength < 0!!  tdex = %d\n",jij);
      fflush(stderr);
      exit(32);
    }
    // if (jij==tdex-1) exit(32);  // An exit for testing purposes.

    // Cartesian coordinates in SC-3 of the VMP.
    xx = nrpt[0] + ttmp * unit[0];
    yy = nrpt[1] + ttmp * unit[1];
    zz = nrpt[2] + ttmp * unit[2];
    // Heliocentric height of the VMP.
    r = sqrt(xx*xx + yy*yy + zz*zz);

    // Now compute "ardex": the column-index of the voxel.
    // The calculation is done for the three geometries.

#ifdef CYLINDRICAL

    rr = sqrt(xx * xx + yy * yy);
    phiphi = atan2(yy, xx);
    if (phiphi < 0.0)
      phiphi += 2.0 * M_PI;

    index[0] = floor((rr / rmax) * (double) NRAD);
    index[1] = floor((phiphi / 2.0 / M_PI) * (double) NPHI);
    index[2] = floor((zz + rmax) / deltagrid);

    ardex = index[2]*NRAD*NPHI + index[1]*NRAD + index[0];

#elif defined HOLLOW_SPHERE

    /* r index */
   //index[0] = floor((double) NRAD *(r - (double) RMIN)/(rmax - (double) RMIN)) ;
   // previous line now changed to:
    index[0] = rad_bin_number(r);
    if (index[0] == NRAD)
      index[0] = NRAD - 1;   /*these index corrections are due to finite precision */

    // theta (polar) index,
    // rr: polar angle, taking rr=0 in South Pole, rr=Pi in North Pole.
    rr = atan( zz / sqrt(xx*xx + yy*yy) ) + M_PI/2.;
    index[1] = floor( rr*((double) NTHETA)/M_PI );
    if (index[1] == NTHETA)
      index[1] = NTHETA - 1;

    /* phi (azimuthal) index */
    // phiphi: longitude, making sure is in the range [0,2Pi]
    phiphi = atan2(yy, xx);
    if (phiphi < 0.0)
      phiphi += 2.0 * M_PI;
    index[2] = floor((phiphi / 2.0 / M_PI) * (double) NPHI);
    if (index[2] == NPHI)
      index[2] = NPHI - 1;

    // In the following formula for "ardex", voxels are numbered
    // starting at [irad,ith,iph]=[0,0,0] and then ordered:
    // first by rad, then by theta, finally by phi.
    ardex = index[2]*NRAD*NTHETA + index[1]*NRAD + index[0];

#elif defined CARTESIAN

    index[0] = floor((xx + rmax) / deltagrid);
    index[1] = floor((yy + rmax) / deltagrid);
    index[2] = floor((zz + rmax) / deltagrid);

    ardex = index[2]*NCELLS*NCELLS + index[1]*NCELLS + index[0];

#endif

    if (ardex > NBINS - 1){
      fprintf(stderr,"ardex (%d) exceeds max. voxel number (%d)!\n", ardex,NBINS);
      fprintf(stderr,"voxel index = (%d,%d,%d), location = (%g,%g,%g)\n",
	      index[0],index[1],index[2],xx,yy,zz);
      fflush(stderr);
      exit(32);
    }

    // At this point we have the VMP "time" (ttmp), the VAL (arclength), and ardex.

    vdhC = 0.;
    vdhD = 0.;
#ifdef THOMSON     /* Thomson scattering */
    rtmp = impact * impact / r / r;
    sgam = 1.0 / r;
    cgam = sqrt(1.0 - sgam * sgam);
    vdhA = cgam * sgam * sgam;
    vdhB = 3. * sgam * sgam - 1.0;
    vdhB +=  /* note: 4-3*cgam^2 = 1 + 3*sgam^2  */
      (cgam*cgam/sgam)*(4. - 3.*cgam*cgam)*log((1.+sgam)/cgam);
    vdhB /= 8.0;
    if (totalB == 0){
      Arow_long[ardex] += (float) ((1. - QLIMB) * vdhA + QLIMB * vdhB) *
                          rtmp * arclength * CONST * RSUN * 1.0e5;
    } else if (totalB == 1) {
      vdhC = 4./3. - cgam*(1 + cgam*cgam/3.);
      vdhD = 5. + sgam*sgam - (cgam*cgam/sgam)*(5. - sgam*sgam)*log((1.+sgam)/cgam);
      vdhD /= 8.0;
      Arow_long[ardex] += (float) ( (1. - QLIMB)*(2.*vdhC - vdhA*rtmp)
  				     +     QLIMB*(2.*vdhD - vdhB*rtmp) )*
                     arclength*CONST*RSUN*1.0e5;
    } else {
      fprintf(stderr,"buildrow.c: totalB = %d!\n",totalB);
      exit(1);
    }
#elif defined RADON
    /* Unweighted ray transform.  The path length is in units of Rs. */
    Arow_long[ardex] += (float) arclength ;
#elif defined PROJTEST
    /* make up your own weighting function!  What a country! */
#elif defined THOMSON_PT
    /*  Thomson scatter from a point source */
    rtmp = impact * impact / r / r;
    vdhA = 1.0 / r / r;
    Arow_long[ardex] +=
      (float) vdhA *rtmp * arclength * CONST * RSUN * 1.0e5;
#endif

#ifdef COMPARISON_CALCULATION /* this is for compare.c */
    pBcalc[ll][kk] += x[ardex]*Arow_long[ardex];
#endif

#ifdef RAYDIAGNOSE
    fprintf(stderr,"(%d,%d,%d, %d)",index[0],index[1],index[2],ardex);
    fflush(stderr);
#ifdef GEOMTEST_OUT /*see geomtest.c */
    fwrite(&ttmp, sizeof(double), 1, fid_geomtest);
    ttmp = (double) Arow_long[ardex]/(arclength*RSUN*1.e5);
    fwrite(&ttmp, sizeof(double), 1, fid_geomtest);
#endif
#endif


  } /*tdex loop */

salida:/* exit point for LOS's that miss the compution grid */

  if (ontarget == 0)
    hasdata = 0;

#ifdef RAYDIAGNOSE
  if (ontarget == 0){
    fprintf(stderr,"EXITING BUILDROW (ontarget = 0) \n");
    fflush(stderr);
  } else {
    fprintf(stderr,"\nEXITING BUILDROW (ontarget = 1) \n");
    fflush(stderr);
  }
#endif

  fflush(stderr);

}/************ END BUILDROW ***********************/
