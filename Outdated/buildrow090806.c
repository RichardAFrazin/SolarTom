/* WARNING!!!
 *
 * Any change made in this file must also be made in compare.c!
 *
 */

/*see geomtest.c uses RAYDIAGNOSE definition */ 

/* #define RAYDIAGNOSE */

/***********  BEGIN BUILDROW ***********************/
{
  double nrpt[3], g1[3], unit[3], los1[3], los2[3];
  double t1, t2, arclength, xx, yy, zz, impact, r, vdhA, vdhB;
  double deltagrid, grideps, rayeps;
  double junk[6], t[NBINS], rtmp, ttmp, gam, sgam, cgam, ptmp;
  int binbin[6], jij, tdex, index[3], ardex, ontarget;

  rayeps = 1.e-12;
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
  fprintf(stderr,"ENTERING BUILDROW: rho1 = %g, eta1 = %g\n",rho1, eta1);
  fflush(stderr);
#endif

  /* unit is the LOS unit vector
   * nrpt is the nearest point vector */

  /* this is correct if eta1 is the usual solar PA (in radians) */
  Rx = rotx(eta1);

  r3tmp[0] = sin(rho1) * sin(rho1) / cos(rho1);
  r3tmp[1] = 0.0;
  r3tmp[2] = sin(rho1) ;
  rotvmul(nrpt, Rx, r3tmp);
  r3scalmul(nrpt, dsun);

  g1[0] = -cos(rho1);
  g1[1] = 0.0;
  g1[2] = sin(rho1);
  rotvmul(unit, Rx, g1);

  free(Rx);

  rotvmul(r3tmp, &R23, nrpt);
  r3eq(nrpt, r3tmp);
  rotvmul(r3tmp, &R23, unit);
  r3eq(unit, r3tmp);

  impact = sqrt(r3dot(nrpt, nrpt));

  /* calculate t1,t2, the "times" where ray enters and leaves computation region
   * los1 and los2 mark where the LOS enter and leave the computation area  */

  /* this is for the 3 cartesian coordinates
      the required components of junk[] are overwritten 
      for the other geometries  */
  for (jij = 0; jij < 6; jij += 2) {
     /*these are the (signed) distances from nrpt where the LOS crosses the max and min
	   values of the computation box for each of the 3 cartesian coordinates */
    junk[jij] = (-rmax - nrpt[jij / 2]) / (unit[jij / 2] + grideps);
    junk[jij + 1] = (rmax - nrpt[jij / 2]) / (unit[jij / 2] + grideps);
  }

  /* If the LOS misses the computation grid, exit buildrow.c .
   * If any of the endpoints specified by the junk vector
   *    hits the edge of the grid, set ontarget = 1 */


#ifdef CARTESIAN
  ontarget = 0;
  for (jij = 0; jij < 6; jij++){
    g1[0] = nrpt[0] + junk[jij]*unit[0];
    g1[1] = nrpt[1] + junk[jij]*unit[1];
    g1[2] = nrpt[2] + junk[jij]*unit[2];

    if ( (fabs(g1[0]) < rmax + grideps) &&
         (fabs(g1[1]) < rmax + grideps) &&
         (fabs(g1[2]) < rmax + grideps) ) {
            ontarget = 1;
#ifdef  RAYDIAGNOSE
      fprintf(stderr,"intersection with computation cube face %d\n",i);
#endif 
    }
  }
  if (ontarget == 0)
    goto salida;

#endif

#ifdef CYLINDRICAL

  vdhA = (unit[0] * unit[0] + unit[1] * unit[1]);
  gam = 2.0 * (unit[0] * nrpt[0] + unit[1] * nrpt[1]);
  rtmp = nrpt[0] * nrpt[0] + nrpt[1] * nrpt[1] - rmax * rmax;
     /* find distances where it enters and leaves the 
	       infinite cylinder with radius rmax */

  junk[0] =
    (-gam - sqrt(gam * gam - 4. * vdhA * rtmp)) / (2. * vdhA + rayeps);
  junk[1] =
    (-gam + sqrt(gam * gam - 4. * vdhA * rtmp)) / (2. * vdhA + rayeps);
  junk[2] = -10.*rmax;	/* these will not contribute */
  junk[3] =  10.*rmax;	

 /* if any of the endpoints specified by the junk vector
     hits the cylinder boundaries, set ontarget = 1 */

#ifdef  RAYDIAGNOSE
        fprintf(stderr,"grideps = %1.10g, rmax+gridesp = %1.10g\n",grideps,rmax + grideps );
        fflush(stderr);
#endif 


  ontarget = 0;
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
#ifdef  RAYDIAGNOSE
            fprintf(stderr,"intersection with computation cylinder face %d\n",jij);
            fflush(stderr);
#endif 
    }
  }  
  if ( ontarget == 0)
      goto salida;

#elif defined HOLLOW_SPHERE

  ontarget = 1;
  if (impact > rmax + grideps ){
    /* Is the LOS outside of computation sphere?  If so, just
         treat it has having no data (see build_subA.c)*/
    ontarget = 0;    
    goto salida;
  }

    /* the LOS enters the sphere at the point nrpt + unit*gam,
       where gam = sqrt(rmax^2 - nrpt'*nrpt), 
	   unit points AWAY from the observer, TOWARDS the Sun
	   (recall, in coordinate system 2, the x-axis points
	      TOWARD the observer  
     *
     * los first enters the sphere at this (signed) distance from nprt */

  junk[0] = - sqrt(rmax*rmax - impact*impact);
  /*does the LOS hit the inner sphere (hollow part)? */
  if (impact <= ((double) RMIN)){
	 junk[1] = - sqrt(((double) RMIN)*((double) RMIN) - impact*impact);
   } else {
	 junk[1] =  sqrt(rmax*rmax - impact*impact);
   }

   t1 = junk[0] + grideps;
   t2 = junk[1] - grideps;

#endif


#if defined CARTESIAN || defined CYLINDRICAL
    /* t2 corresponds to the minimum of the positive times
	   t1                    maximum        negative         */
  t1 = -10.0 * rmax;
  t2 = 10.0 *  rmax;
  for (jij = 0; jij < 6; jij++) {
    if (junk[jij] >= grideps) {
      if (junk[jij] < t2) {
        t2 = junk[jij];
      }
    }
    if (junk[jij] < grideps) {
      if (junk[jij] > t1) {
        t1 = junk[jij];
      }
    }
  }
  t1 += grideps;
  t2 -= grideps;
#endif

  for (jij = 0; jij < 3; jij++) {
    los1[jij] = nrpt[jij] + t1 * unit[jij];
    los2[jij] = nrpt[jij] + t2 * unit[jij];
  }

  /* put the bin number of LOS endpoints into binbin array -
       see CARTESIAN example for ordering */

#ifdef CYLINDRICAL

  rtmp = sqrt(los1[0] * los1[0] + los1[1] * los1[1]);
  binbin[0] = floor(rtmp * ((double) NRAD) / rmax);
  rtmp = sqrt(los2[0] * los2[0] + los2[1] * los2[1]);
  binbin[1] = floor(rtmp * ((double) NRAD) / rmax);

  ptmp = atan2(los1[1], los1[0]);
  if (ptmp < 0.0)
    ptmp += 2.0 * M_PI;
  binbin[2] = floor(ptmp * ((double) NPHI) / 2.0 / M_PI);
#ifdef RAYDIAGNOSE
  fprintf(stderr,"entry phi = %g deg, ",ptmp*180./M_PI);
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
  fprintf(stderr,"exit phi = %g deg |(entry phi) - (exit phi)| = %g, ==> wrap = %d\n", 
    rtmp*180./M_PI,fabs(rtmp - ptmp)*180./M_PI,wrap);
#endif

  binbin[4] = floor((los1[2] + rmax) / deltagrid);
  binbin[5] = floor((los2[2] + rmax) / deltagrid);

#elif defined CARTESIAN

  binbin[0] = floor((los1[0] + rmax) / deltagrid);
  binbin[1] = floor((los2[0] + rmax) / deltagrid);
  binbin[2] = floor((los1[1] + rmax) / deltagrid);
  binbin[3] = floor((los2[1] + rmax) / deltagrid);
  binbin[4] = floor((los1[2] + rmax) / deltagrid);
  binbin[5] = floor((los2[2] + rmax) / deltagrid);

#elif defined HOLLOW_SPHERE
      /* r, theta (polar angle), phi (azimuthal angle) coordinate order */
  binbin[0] = NRAD-1;
  binbin[1] = NRAD-1;
  if (impact <= (double) RMIN){
     binbin[1] = 0;
  }
    
  rtmp = atan( los1[2] / sqrt( los1[0]*los1[0] + los1[1]*los1[1] + rayeps) );
  binbin[2] =  floor( (rtmp + M_PI/2.0)*((double) NTHETA)/ M_PI );
  rtmp = atan( los2[2] / sqrt( los2[0]*los2[0] + los2[1]*los2[1] + rayeps) );
  binbin[3] =  floor( (rtmp + M_PI/2.0)*((double) NTHETA)/ M_PI );
  
  rtmp = atan2(los1[1], los1[0]);
  if (rtmp < 0.0) rtmp += 2.0 * M_PI;
  binbin[4] = floor(rtmp * ((double) NPHI) / 2.0 / M_PI);
#ifdef RAYDIAGNOSE    
  fprintf(stderr,"entry phi = %g deg, ",rtmp*180./M_PI);
#endif
  ptmp = atan2(los2[1], los2[0]);
  if (ptmp < 0.0) ptmp += 2.0 * M_PI;
  binbin[5] = floor(ptmp * ((double) NPHI) / 2.0 / M_PI);
  wrap = 0;
  if ( fabs(rtmp - ptmp) > M_PI )
        wrap = 1;
#ifdef RAYDIAGNOSE  
  fprintf(stderr,"exit phi = %g deg, |(entry phi) - (exit phi)| = %g ===> wrap = %d\n",
      ptmp*180/M_PI, fabs(rtmp - ptmp)*180./M_PI, wrap);
#endif

#endif

  /* find the "times" of bin crossings */

  t[0] = t1;
  t[1] = t2;
  tdex = 2;

#ifdef RAYDIAGNOSE  
  fprintf(stderr,"entry distance t1 = %g, exit distance t2 = %g,\n",t1,t2);
  fprintf(stderr,"nrpt = %g %g %g\n",nrpt[0],nrpt[1],nrpt[2]);
  fprintf(stderr,"unit = %g %g %g\n",unit[0],unit[1],unit[2]);
  fprintf(stderr,"point of entry %g, %g, %g\n",nrpt[0]+ t1*unit[0], nrpt[1]+t1*unit[1], nrpt[2]+t1*unit[2]);
  fprintf(stderr,"point of exit  %g, %g, %g\n",nrpt[0]+ t2*unit[0], nrpt[1]+t2*unit[1], nrpt[2]+t2*unit[2]);
  fprintf(stderr,"entry/exit bins (binbin) = %d %d %d %d %d %d\n",binbin[0],binbin[1],binbin[2],binbin[3],binbin[4],binbin[5]);
#endif


#ifdef CYLINDRICAL

  /* binrmin = bin of minimum radius */
  ttmp = -(unit[0] * nrpt[0] + unit[1] * nrpt[1]) /
         (unit[0] * unit[0] + unit[1] * unit[1] + rayeps);
  rtmp = sqrt((nrpt[0] + ttmp * unit[0]) * (nrpt[0] + ttmp * unit[0]) +
              (nrpt[1] + ttmp * unit[1]) * (nrpt[1] + ttmp * unit[1]));
  binrmin = floor(rtmp * ((double) NRAD) / rmax + rayeps);

  rtmp = nrpt[2] + ttmp * unit[2];

  /*radial bin crossings */
  vdhA = (unit[0] * unit[0] + unit[1] * unit[1]);
  gam = 2.0 * (unit[0] * nrpt[0] + unit[1] * nrpt[1]);
  ptmp = nrpt[0] * nrpt[0] + nrpt[1] * nrpt[1];

#ifdef RAYDIAGNOSE
  fprintf(stderr,"RADIAL Bins: binrmin= %d, bin crossings: ",binrmin);
  fflush(stderr);
#endif

  for (jij = MAX(binbin[0], binbin[1]); jij >= binrmin; jij--) {
    rtmp = rmax * ((double) (jij + 1) / (double) NRAD);
    rtmp = ptmp - rtmp * rtmp;
    ttmp =
      (-gam - sqrt(gam * gam - 4. * vdhA * rtmp)) / (2. * vdhA +
          rayeps);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d,%g)",jij,ttmp);
      fflush(stderr);
#endif
    } 
    ttmp =
      (-gam + sqrt(gam * gam - 4. * vdhA * rtmp)) / (2. * vdhA +
          rayeps);
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

  /*  does the LOS go through the 1/2 plane y = 0, x > 0
   *  within the computation boundary?   see above */
 
  /* old wrapping calculation.  Seems to be  
       OK but the method above is simpler
 
      wrap = 0;
      ttmp = -nrpt[1] / (rayeps + unit[1]);
      rtmp = nrpt[0] + ttmp * unit[0];
      gam = nrpt[2] + ttmp * unit[2];
      if (((rtmp > 0) && (rtmp <= rmax))
         || ((gam >= -rmax) && (gam <= rmax)))
      wrap = 1;
  */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nPHI Bins: bin crossings: ");
  fflush(stderr);
#endif

  if (wrap == 0) {
    for (jij = MIN(binbin[2], binbin[3]) + 1;
         jij <= MAX(binbin[2], binbin[3]); jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] +
                                      rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 0: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij * 360.0/ (double) NPHI));
        exit(32);
      }
    }
  } else {
    for (jij = MAX(binbin[2], binbin[3]) + 1; jij < NPHI; jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] +
                                      rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij * 360.0/ (double) NPHI));
        exit(32);
      }
    }
    for (jij = 0; jij <= MIN(binbin[2], binbin[3]); jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] +
                                      rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,ttmp);
        fflush(stderr);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
        fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij * 360.0/ (double) NPHI));
        exit(32);
      }
    }
  }

  /* z bin crossings */

#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nZ Bins: bin crossings: ");
  fflush(stderr);
#endif
  
  for (jij = MIN(binbin[4], binbin[5]); jij < MAX(binbin[4], binbin[5]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[2]) / (unit[2] + rayeps);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
#ifdef RAYDIAGNOSE
      fprintf(stderr,"(%d,%g)",jij,ttmp);
      fflush(stderr);
#endif
    } else {
      fprintf(stderr,"ttmp = %g, ptmp = %g, jij = %d, z = %g \n",ttmp,ptmp, jij, (jij+1)*deltagrid - rmax);
      exit(32);
    }
  }
  
#elif defined HOLLOW_SPHERE
 /*radial bin crossings */
	
  if (impact <= ((double) RMIN)){
	 binrmin = 0;  /* binrmin = bin of minimum radius */
   } 
   else {
	 binrmin = floor( (impact - (double) RMIN) * ((double) NRAD) / (rmax - (double) RMIN)) ;
   }
    /* vdhA is \Delta r */
  vdhA = (rmax - ((double) RMIN)) / ((double) NRAD);
  for (jij = NRAD - 1; jij >= binrmin; jij--) {
    rtmp = ((double) RMIN) + (jij + 1)*vdhA;
    ttmp = sqrt(rtmp*rtmp - impact*impact) ;
	t[tdex] =   ttmp;
	tdex++;
	if ( impact > (double) RMIN ) {
  	  t[tdex] = - ttmp;
	  tdex++;
	}
  }
  /* polar angle bin crossings */

  for (jij = MIN(binbin[2],binbin[3]); jij <= MAX(binbin[2],binbin[3]); jij++) {
	 gam = tan( jij*M_PI/((double) NTHETA) - M_PI/2.0 );
	 gam = gam*gam;
	 vdhA = gam*(unit[0]*unit[0] + unit[1]*unit[1]) - unit[2]*unit[2];
	 vdhB = 2.*( gam*(unit[0] + unit[1]) - unit[2]);
	 sgam = gam*(nrpt[0]*nrpt[0] + nrpt[1]*nrpt[1]) - nrpt[2]*nrpt[2];
	 
	 ttmp =  (- vdhB + sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA + rayeps);
	 if ((ttmp > t1) && (ttmp < t2)) {
	    t[tdex] = ttmp;
		tdex++;
	 }
     ttmp =  (- vdhB - sqrt(vdhB*vdhB - 4.*vdhA*sgam))/(2.*vdhA + rayeps);
	 if ((ttmp > t1) && (ttmp < t2)) {
	    t[tdex] = ttmp;
		tdex++;
	 }
  }
  /* azimuthal bin crossings */

  /*  Does the LOS go through the 1/2 plane y = 0, x > 0 within the
        computation boundary?  If so, wrap = 1 (see above).  
   */
  


  if (wrap == 0) {
    for (jij = MIN(binbin[4], binbin[5]) + 1;
         jij <= MAX(binbin[4], binbin[5]); jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] + rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE   
        fprintf(stderr,"jij= %d, ttmp= %g  ",jij, ttmp);
#endif
      } else {
        fprintf(stderr, "wrap = 0: out of bounds!!\n");
	fprintf(stderr, "ttmp = %g, ptmp = %g, jij = %d, phi = %g deg\n",ttmp,ptmp, jij,(jij * 360.0/ (double) NPHI));
	exit(32);
      }
    }
  } else {
    for (jij = MAX(binbin[4], binbin[5])+1; jij < NPHI; jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] +
                                      rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE 
        fprintf(stderr,"jij= %d, ttmp= %g  ",jij, ttmp);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n");
	fprintf(stderr, "ttmp = %g, jij = %d\n",ttmp,jij);
	exit(32);

      }
    } /* jij loop */
    for (jij = 0; jij <= MIN(binbin[4], binbin[5]); jij++) {
      ptmp = tan(jij * 2.0 * M_PI / (double) NPHI);
      ttmp =
        (nrpt[1] - nrpt[0] * ptmp) / (unit[0] * ptmp - unit[1] +
                                      rayeps);
      if ((ttmp > t1) && (ttmp < t2)) {
        t[tdex] = ttmp;
        tdex++;
#ifdef RAYDIAGNOSE 
        fprintf(stderr,"jij= %d, ttmp= %g  ",jij, ttmp);
#endif
      } else {
        fprintf(stderr, "wrap = 1: out of bounds!!\n"); 
        fprintf(stderr,"jij= %d, ttmp= %g  ",jij, ttmp);
	exit(32); 
      }
    } /* jij loop */
  }

   
#elif defined (CARTESIAN)

  /* x bin crossings */

  for (jij = MIN(binbin[0], binbin[1]); jij < MAX(binbin[0], binbin[1]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[0]) / (unit[0] + rayeps);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }

  /* y bin crossings */

  for (jij = MIN(binbin[2], binbin[3]); jij < MAX(binbin[2], binbin[3]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[1]) / (unit[1] + rayeps);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }

  /* z bin crossings */

  for (jij = MIN(binbin[4], binbin[5]); jij < MAX(binbin[4], binbin[5]);
       jij++) {
    ptmp = (jij + 1) * deltagrid - rmax;
    ttmp = (ptmp - nrpt[2]) / (unit[2] + rayeps);
    if ((ttmp > t1) && (ttmp < t2)) {
      t[tdex] = ttmp;
      tdex++;
    }
  }

#endif


  /* sort the "times" */
  qsort((void *) t, tdex, sizeof(double), (void *) &doublecompare);


#ifdef RAYDIAGNOSE
  fprintf(stderr,"\nVoxel Numbers: ");
  fflush(stderr);
#endif

  for (jij = 1; jij < tdex; jij++) {

    /* calculate the matrix elements */

    ttmp = 0.5 * (t[jij] + t[jij - 1]);
    arclength = t[jij] - t[jij - 1];

    if (arclength < 0.0) {
      fprintf(stderr, "BUILDROW: arclength < 0!!  tdex = %d\n",jij);
      fflush(stderr);
      exit(32);
    }


    xx = nrpt[0] + ttmp * unit[0];
    yy = nrpt[1] + ttmp * unit[1];
    zz = nrpt[2] + ttmp * unit[2];
    r = sqrt(xx*xx + yy*yy + zz*zz);

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

    index[0] = floor((r / rmax) * (double) NRAD);

    rr = atan( zz / ( sqrt(xx*xx + yy*yy) + grideps) ) + M_PI/2.;
    index[1] = floor( rr*(double (NTHETA))/M_PI );

    phiphi = atan2(yy, xx);
    if (phiphi < 0.0)
      phiphi += 2.0 * M_PI;
    index[2] = floor((phiphi / 2.0 / M_PI) * (double) NPHI);

    ardex = index[2]*NRAD*NTHETA + index[1]*NRAD + index[0];

#elif defined CARTESIAN

    index[0] = floor((xx + rmax) / deltagrid);
    index[1] = floor((yy + rmax) / deltagrid);
    index[2] = floor((zz + rmax) / deltagrid);
    ardex = index[2]*NCELLS*NCELLS + index[1]*NCELLS + index[0];

#endif

    if (ardex > NBINS - 1){
      fprintf(stderr,"ardex (%d) exceeds max. voxel number (%d)!\n",ardex,NBINS);
      fprintf(stderr,"location = (%g,%g,%g)\n",xx,yy,zz);
      fflush(stderr);
      exit(32);
    }   

#ifdef THOMSON
    /* Thomson scattering */
    rtmp = impact * impact / r / r;
    sgam = 1.0 / r;
    cgam = sqrt(1.0 - sgam * sgam);
    vdhA = cgam * sgam * sgam;
    vdhB = 3. * sgam * sgam - 1.0;
    vdhB +=
      (cgam * cgam / sgam) * (4. -
                              3. * cgam * cgam) * log((1. +
                                                       sgam) / cgam);
    vdhB /= 8.0;
    Arow_long[ardex] += (float) ((1. - QLIMB) * vdhA + QLIMB * vdhB) *
                        rtmp * arclength * CONST * RSUN * 1.0e5;

#elif defined RADON
    /* do the straight Radon transform with the same scaling */
    Arow_long[ardex] += (float) arclength *CONST * RSUN * 1.0e5;
#elif defined PROJTEST
    /* make up your own weighting function!  What a country! */
#elif defined THOMSON_PT
    /*  Thomson scatter from a point source */
    rtmp = impact * impact / r / r;
    vdhA = 1.0 / r / r;
    Arow_long[ardex] +=
      (float) vdhA *rtmp * arclength * CONST * RSUN * 1.0e5;
#endif

#ifdef RAYDIAGNOSE
    fprintf(stderr," %d",ardex);
    fflush(stderr);
#endif


  } /*tdex loop */

salida: ;  /* exit point for LOS's that miss the compution grid */

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
