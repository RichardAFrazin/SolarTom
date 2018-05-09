/* HOLLOW_SPHERE: (r, theta, phi)
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
  double nrpt[3], g1[3], unit[3], los1[3], los2[3];
  double t1, t2, t3, arclength, xx, yy, zz, impact, r, vdhA, vdhB, vdhC, vdhD;
  double junk[2], t[NBINS], rtmp, ttmp, gam, sgam, cgam, ptmp;
  int binbin[6], jij, tdex, index[3], ardex, ontarget;
  double abstrmin, abstrmax, bin_bdy[2]; 
  int index0, case_num;
  char case_str[]="xx";
  int wrap, binrmin;
  double rr, phiphi;

#ifdef RAYDIAGNOSE
  fprintf(stderr,"ENTERING BUILDROW: rho1 = %1.12g, eta1 = %1.12g\n",rho1, eta1);
  fflush(stderr);
#endif

  junk[0] = 0.; junk[1] = 0.;

  /* unit is the LOS unit vector
   * nrpt is the nearest point vector
   *
   * unit points AWAY from the observer, TOWARDS the Sun
   *   (recall, in coordinate system 2, the x-axis points
   *   TOWARD the observer) */

  // Note in CASE-3-D Unit MAY point AWAY from the observer and AWAY the Sun,
  // I will ASSUME this is already possible in the current code
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

  // LOS impact parameter = Norm(NRPT)
  impact = sqrt(r3dot(nrpt, nrpt)); // [Rsun] units

  // Compute t3, the SIGNED "time" of the spacecraft position alog the LOS.
  t3 = sqrt(dsun*dsun - impact*impact);
  if (r3dot(unit,sun_ob3) < 0)
    t3 *= -1.;

  // If the LOS misses the computation grid set ontarget = 0,
  // treating it as if it has no data, and exit builrow.c
  // Note we allow now two possible situations for this.
  ontarget = 1;
  if ( (impact > RMAX) || (dsun < RMIN && impact <= 1. && t3 < 0.) ){
    ontarget = 0;
    goto salida;
  }

  // Figure out case_str for spacecraft position.
  case_num = 0;
  if (dsun > ((double) RMAX)){
    case_str[0] = '1';
    if (impact >= RMIN){
      case_str[1] = 'A';
      case_num    = 11 ;
    } else if (impact < RMIN && impact > 1.){
        case_str[1] = 'B';
        case_num    = 12 ;
    } else if (impact < 1.){
        case_str[1] = 'C';
        case_num    = 13 ;
    } else {
        fprintf(stderr, "buildrow: unrecognized case 1 of spacecraft position.");
        fprintf(stderr, "RMIN = %g, RMAX = %g, dsun = %g, impact = %g, t3 = %g",RMIN,RMAX,dsun,impact,t3);
        exit(0);
    }      
  } else if ((dsun <= (double) RMAX) && (dsun >= (double) RMIN)){
    case_str[0] = '2';
    if (impact >= RMIN && t3 <= 0.){
        case_str[1] = 'A';
	case_num    = 21 ;
    } else if (impact >= RMIN && t3 > 0.) {
        case_str[1] = 'B';
	case_num    = 22 ;
    } else if (impact < RMIN && impact > 1. && t3 < 0.) {
        case_str[1] = 'C';
	case_num    = 23 ;
    } else if (impact < RMIN && impact > 1. && t3 > 0.) {
        case_str[1] = 'D';
	case_num    = 24 ;
    } else if (impact < 1. && t3 < 0.) { 
        case_str[1] = 'E';
	case_num    = 25 ;
    } else if (impact < 1. && t3 > 0.) { 
        case_str[1] = 'F';
	case_num    = 26 ;
    } else {
        fprintf(stderr, "buildrow: unrecognized case 2 of spacecraft position.");
        fprintf(stderr, "RMIN = %g, RMAX = %g, dsun = %g, impact = %g, t3 = %g",RMIN,RMAX,dsun,impact,t3);
        exit(0);
    } 
  } else { // if (dsun < RMIN)
    case_str[0] = '3';
    if (impact > 1. && t3 <= 0.){ 
        case_str[1] = 'A';
	case_num    = 31 ;
    } else if (impact >  1. && t3 > 0.) {
        case_str[1] = 'B';
	case_num    = 32 ;
    } else if (impact <= 1. && t3 < 0.) {
        case_str[1] = 'C';
	case_num    = 33 ;
	fprintf(stderr, "buildrow: case_str valued 3C was allowed.");
	fprintf(stderr, "RMIN = %g, RMAX = %g, dsun = %g, impact = %g, t3 = %g",RMIN,RMAX,dsun,impact,t3);
	exit(0);
    } else if (impact <= 1. && t3 > 0.) {
        case_str[1] = 'D';
	case_num    = 34 ;
    } else {
        fprintf(stderr, "buildrow: unrecognized case 3 of spacecraft position.");
        fprintf(stderr, "RMIN = %g, RMAX = %g, dsun = %g, impact = %g, t3 = %g",RMIN,RMAX,dsun,impact,t3);
        exit(0);
    }
  }

  //------------------------------------------------------
  /* Calculate t1,t2, the "times" where ray enters and leaves computation region
   * los1 and los2 mark where the LOS enters and leaves the computation area
   *
   * junk[] are the (signed) distances from nrpt where the LOS crosses
   *    the max and min values of the computation ball for each of the 3
   *    coordinates
   *
   * If any of the endpoints specified by the junk vector
   *    hits the edge of the grid, set ontarget = 1
   */

  // Un-signed crossing times at radii RMIN and RMAX.
  abstrmin = sqrt(((double) RMIN)*((double) RMIN) - impact*impact);
  abstrmax = sqrt(((double) RMAX)*((double) RMAX) - impact*impact);

  // For all possible on-target cases assign:
  // t1 and t1, binbin[0] and binbin[1]
  switch (case_num) {
  case 11://"1A":
  case 12://"1B":
  case 21://"2A":
  case 22://"2B":
  case 23://"2C":
  case 24://"2D":
    t1 = -abstrmax;
    t2 =  abstrmax;
    binbin[0] = NRAD-1;
    binbin[1] = NRAD-1;
    break;
  case 13://"1C":
  case 25://"2E":
    t1 = -abstrmax;
    t2 = -abstrmin;
    binbin[0] = NRAD-1;
    binbin[1] = 0;
    break;
  case 26://"2F":
  case 31://"3A":
  case 32://"3B":
  case 34://"3D":
    t1 =  abstrmin;
    t2 =  abstrmax;
    binbin[0] = 0;
    binbin[1] = NRAD-1;
    break;
  default:
    fprintf(stderr, "buildrow: case_str = %s",case_str);
    exit(0);
  }
  
// LOS endpoints position vector coordinates in CS-3.
// This is valid for all geometries.
   for (jij = 0; jij < 3; jij++)
   {
    los1[jij] = nrpt[jij] + t1*unit[jij];
    los2[jij] = nrpt[jij] + t2*unit[jij];
   }

  // Put the bin number of LOS endpoints into binbin array in
  // r, theta (polar angle), phi (azimuthal angle) coordinate order.
  // Components 0 and 1 of binbin already assigned above, do components 2,3,4,5 */

  // Compute Latitude [rad] of vector los1
  rtmp = atan( los1[2] / sqrt( los1[0]*los1[0] + los1[1]*los1[1]) );
  // Compute Theta bin of vector los1, taking theta=0 as the South pole,  theta=Pi as the North pole.
  binbin[2] =  floor( (rtmp + M_PI/2.)*((double) NTHETA)/ M_PI );  // assumes uniform theta grid
  
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry theta = %g deg, ",rtmp*180./M_PI);
#endif
  // Same Latitude/Theta computations for vector los2:
  rtmp = atan( los2[2] / sqrt( los2[0]*los2[0] + los2[1]*los2[1]) );
  binbin[3] =  floor( (rtmp + M_PI/2.)*((double) NTHETA)/ M_PI );
#ifdef RAYDIAGNOSE
  fprintf(stderr,"exit theta = %g deg\n",rtmp*180/M_PI, wrap); // Here wrap has not been assigned a value yet. 
#endif
  // Compute Longitude [rad] of vector los1
  rtmp = atan2(los1[1], los1[0]);
  if (rtmp < 0.)
    rtmp += 2.*M_PI;
  // Compute Phi bin of vector los1
  binbin[4] = floor(rtmp * ((double) NPHI)/2./M_PI);
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry phi = %g deg, ",rtmp*180./M_PI);
#endif
  // Same Longitude/Phi computations for vector los2:
  ptmp = atan2(los2[1], los2[0]);
  if (ptmp < 0.)
    ptmp += 2.*M_PI;
  binbin[5] = floor(ptmp * ((double) NPHI)/2./M_PI);
  wrap = 0;
  if ( fabs(rtmp - ptmp) > M_PI )  // check for wrapping (crossing phi = 0) 
    wrap = 1; 
#ifdef RAYDIAGNOSE
  fprintf(stderr,"exit phi = %g deg, ==> wrap = %d\n",ptmp*180/M_PI, wrap);
#endif

  /* find the "times" of bin crossings */

  // First two elements are entry and exit "times"
  t[0] = t1;
  t[1] = t2;
  tdex = 2;

  // Third element is the Spacecraft "time"
  t[2] = t3;
  tdex++;

#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry distance t1 = %g, exit distance t2 = %g,\n",t1,t2);
  fprintf(stderr,"nrpt = %1.10g %1.10g %1.10g\n",nrpt[0],nrpt[1],nrpt[2]);
  fprintf(stderr,"unit = %1.10g %1.10g %1.10g\n",unit[0],unit[1],unit[2]);
  fprintf(stderr,"point of entry %1.10g, %1.10g, %1.10g\n",nrpt[0]+ t1*unit[0], nrpt[1]+t1*unit[1], nrpt[2]+t1*unit[2]);
  fprintf(stderr,"point of exit  %1.10g, %1.10g, %1.10g\n",nrpt[0]+ t2*unit[0], nrpt[1]+t2*unit[1], nrpt[2]+t2*unit[2]);
  fprintf(stderr,"entry/exit bins (binbin) = %d %d %d %d %d %d\n",binbin[0],binbin[1],binbin[2],binbin[3],binbin[4],binbin[5]);
#endif

  /*radial bin crossings */
  if (impact <= ((double) RMIN)){
     binrmin = 0;  /* binrmin = bin of minimum radius */
  } else {
     binrmin = rad_bin_number(impact);
  }
#ifdef RAYDIAGNOSE
  fprintf(stderr,"RADIAL Bins: binrmin= %d, bin crossings: ",binrmin);
  fflush(stderr);
#endif
  
  /* -2 because of the bin numbering and the entry point into the last bin is already marked by t1 */
  // or by t2 (case 2F and all cases 3).
  for (jij = NRAD - 2; jij >= binrmin; jij--) {
    rad_bin_boundaries(jij, bin_bdy);  
    rtmp = *bin_bdy; // outer boundary of cell jij.

    // Compute here the absolute value of the crossing time at rtmp, naming it ttmp.
    // For cases 1, 2, 3, decide when to assign it negative sign or positive sign, or take both when appropriate.
   ttmp = sqrt(rtmp*rtmp - impact*impact) ;

   switch (case_num) {
   case 11:
   case 12:
   case 21:
   case 22:
   case 23:
   case 24:
     t[tdex] = -ttmp;
#ifdef RAYDIAGNOSE
     fprintf(stderr,"(%d,%g)",jij,t[dex]);
     fflush(stderr);
#endif
     tdex++;
     t[tdex] = +ttmp;
#ifdef RAYDIAGNOSE
     fprintf(stderr,"(%d,%g)",jij,t[dex]);
     fflush(stderr);
#endif
     tdex++;
     break;
   case 13:
   case 25:	  
     t[tdex] = -ttmp;
#ifdef RAYDIAGNOSE
     fprintf(stderr,"(%d,%g)",jij,t[dex]);
     fflush(stderr);
#endif
     tdex++;
     break;
   case 26:
   case 31:
   case 32:
   case 34:
     t[tdex] = +ttmp;
#ifdef RAYDIAGNOSE
     fprintf(stderr,"(%d,%g)",jij,t[dex]);
     fflush(stderr);
#endif
     tdex++;
     break;
   default:
     fprintf(stderr, "buildrow: case_str = %s",case_str);
     exit(0);
   }
   
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
    //ttmp =  (- vdhB - sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA);
	ttmp =  GridDivision(- vdhB - sqrt(vdhB*vdhB - 4.*vdhA*sgam), 2.*vdhA);
		 
	if ((ttmp > t1) && (ttmp < t2)) {
	  t[tdex] = ttmp;
	  tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d, %g deg, %g)", jij, (jij+1)*180./((double) NTHETA) - 90., ttmp);
            fflush(stderr);
#endif
    }  

    //ttmp =  (- vdhB + sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA);
    ttmp =  GridDivision(- vdhB + sqrt(vdhB*vdhB - 4.*vdhA*sgam), 2.*vdhA);
    if ((ttmp > t1) && (ttmp < t2)) {
	  t[tdex] = ttmp;
	  tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d, %g deg, %g)",jij,(jij+1)*180./((double) NTHETA) - 90.,ttmp);
      fflush(stderr);
#endif
    }
  }

  /* azimuthal bin crossings */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nPHI Bins: bin crossings: ");
  fflush(stderr);
#endif
  // Here wrap is used in hollow sphere  
  if (wrap == 0) {
    for (jij = MIN(binbin[4], binbin[5]);
         jij < MAX(binbin[4], binbin[5]); jij++) {
      // Let ptmp = TAN(longitude) = r1/r0, where r is the position vector of a cell boundary crossed by the LOS.
      ptmp = tan((jij+1)*2.*M_PI / (double) NPHI);
      // The expression for the crossing time ttmp below comes from solving components 0 and 1 of Eq (8) of FJ2002.
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
    // end of wrap = 0 condition
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
  } // end of wrap = 1 condition


  // At this point, tdex = # of voxels that are threaded by the LOS, +1 (spacecraft location).
 
  // Sort all the bin crossing "times", in ascending order.
  qsort((void *) t, tdex, sizeof(double), (void *) &doublecompare);

  // Find index0 such that t[index0]=t3
  index0 = 0;
  while (t[index0] < t3) index0++;
  //  fprintf(stderr, "t[%d] = %g.  t3 = %g.  t[%d] = %g.\n",index0,t[index0],t3,index0+1,t[index0+1]);
  
  // If case_str[0] NE "2" then index0 = index0 + 1, to exclude spacecraft location.
  // OLD CODE: if (strcmp(case_str,"1") == 0 || strcmp(case_str,"3") == 0) index0++;
  switch (case_num) {
  case 11:
  case 12:
  case 13:
  case 31:
  case 32:
  case 34:
    index0++;
    break;
  default:
    fprintf(stderr, "");
  }
  //  fprintf(stderr, "Case # %s. Set index0 = %d.\n",case_str,index0);
  
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
