
/*---------------code I did for buld_subA.c------------------------*/

for ( i = 0;  i < imsize;  i++) {
    for (jj = 0; jj < imsize; jj++) {
	/* Keep only data within certain radius range  */
#if (defined C2BUILD || defined C3BUILD || defined CORBUILD)
        if (( tan(ARCSECRAD * rho[i][jj]) * dist > INSTR_RMAX * RSUN ) ||
	    ( tan(ARCSECRAD * rho[i][jj]) * dist < INSTR_RMIN * RSUN ))
	pBval[i][jj] = -999.0;
#endif
#if (defined EITBUILD || defined EUVIBUILD || defined AIABUILD)
	if ( (tan(ARCSECRAD * rho[i][jj]) * dist  > INSTR_RMAX * RSUN )
	pBval[i][jj] = -999.0;
#endif
#ifdef RING_REJECT
	if  ((tan(ARCSECRAD * rho[i][jj]) * dist  > INNER_REJECT_RAD * RSUN) &&
	     (tan(ARCSECRAD * rho[i][jj]) * dist  < OUTER_REJECT_RAD * RSUN)) 
	pBval[i][jj] = -999.0;
#endif
#if (defined C2BUILD || defined C3BUILD || defined CORBUILD)
          /* the .79 factor is to convert from units of mean brightness
           *    to 1.e10*(center brightness)           */
if ( abs(pBval[i][jj] + 999.) > QEPS)
  {  /* check for -999 values (missing blocks) */
   #ifdef NRL
   pBval[i][jj] *=	 1.e10 * 0.79;
   #endif
   #ifdef MARSEILLES
   pBval[i][jj] *=  	 0.79; /* Marseilles scaling */
   #endif
  }
#endif
#ifdef DROP_NEG_PB
          if (pBval[i][jj] < 0)
            pBval[i][jj] = -999.0;
#endif
    } /* jj loop over image pixels */
    } /*  i loop over image pixels */  

    /*-------------------------------------------------------------------*/
    
