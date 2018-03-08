
  ttmp = sqrt(rtmp*rtmp - impact*impact) ; // take here the absolute value of the time, and assign proper signs below.
  if (dsun > ((double) RMAX)){ // Cases 1.
      if (impact > 1.0){         // Cases 1A, 1B.
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
       }
     if (impact <= 1.0){        // Case 1C.
       t[tdex] = -ttmp;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,t[dex]);
        fflush(stderr);
#endif
       tdex++;
       }
   } // Cases 1.

   if (dsun >= ((double) RMIN) && dsun <= ((double) RMAX)){ // Cases 2.
      if (impact > 1.0){              // Cases 2A, 2B, 2C, 2D.
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
       }
     if (impact <= 1.0){             // Cases 2E, 2F.
       if (r3dot(unit,sun_ob3) < 0){ // Case 2E.
       t[tdex] = -ttmp;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,t[dex]);
        fflush(stderr);
#endif
       tdex++;
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 2F.
       t[tdex] = +ttmp;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,t[dex]);
        fflush(stderr);
#endif
       tdex++;
       }
       }
   } // Cases 2.

   if (dsun < ((double) RMIN)){      // Cases 3.
      if (impact > 1.0){              // Cases 3A, 3B.
       t[tdex] = +ttmp;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,t[dex]);
        fflush(stderr);
#endif
       tdex++;
       }
     if (impact < 1.0){
       if (r3dot(unit,sun_ob3) < 0){ // Case 3C; set ontarget=1 (LOS hits Sun w/o intersecting grid)
         ontarget = 0;
         goto salida;
       }
       if (r3dot(unit,sun_ob3) > 0){ // Case 3D.
       t[tdex] = +ttmp;
#ifdef RAYDIAGNOSE
        fprintf(stderr,"(%d,%g)",jij,t[dex]);
        fflush(stderr);
#endif
       tdex++;
       }
       }
   } // Cases 3
