/* to compile: 
 * gcc geomtest.c r3misc.o rots.o grids.o -lm -o geomtest
 * note the symbolic link to this file in the directory above
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdlib.h>
#include "headers.h"

#define RAYDIAGNOSE   /* activate print statements in buildrow.c */
#define GEOMTEST_OUT  /* set for buildrow.c to create binary output file
		       * this is meant for comparison to scatter_calc.m 
		       *   read_geomtest.m reads this file.  */

int main( int argc, char **argv)
{
	double sun_ob1[3], sun_ob2[3], spol1[3], sob[3], spol2[3],r3tmp[3];
	double dsun, pang, rho1, eta1, carlong; //deltagrid,
	Rot R12, R23, Rtmp;
	Rot *Rx, *Ry, *Rz;
	double sun_ob3[3];// Extra variables added by Albert, mainly for testing purposes, but also sun_ob3 serves to determine sign of t3, the "time" of the spacecraft.
	int i, hasdata, totalB;
	const double rmax = RMAX;
        float Arow_long[NBINS]; /* for consistency with buildrow.c */
#ifdef GEOMTEST_OUT
	FILE *fid_geomtest;
	char geom_fname[] = "geomtest.dat";/*binary output filename*/
#endif

	totalB = 0; /* for testing scattering calculation */

	//------------These specific values are to be used by Albert-----------------
	sun_ob1[0] = 0.*RSUN;  //so1=[0, 9.04594, 4.26275] Rs
        sun_ob1[1] = 9.0459371*RSUN;
        sun_ob1[2] = 4.2627482*RSUN;
	carlong    = M_PI; /* in radians */
	/* enter the projected radius of the ray in Rsun */
	rho1 =  9.75149;
	/* enter the position angle of ray in solar image */
	eta1 = (270.-0.) * M_PI/180.; 
	//---------------------------------------------------------------------------
	
#ifdef GEOMTEST_OUT
	if ((fid_geomtest = fopen(geom_fname,"wb")) == NULL){
	  fprintf(stderr,"geomtest.c: can't open file: %s\n",geom_fname);
	  exit(1);
	}
#endif
	dsun = 0.0;
	for (i = 0; i < 3; i++) {
		sun_ob1[i] /= RSUN;
		dsun += sun_ob1[i] * sun_ob1[i];
	}
	dsun = sqrt(dsun);

	/* calculate the angle  (NOT IN ORIGINAL CODE) */
	//OLD CODE BY RICH----------------------------------------------------------
	//rho1 = atan( rho1 / dsun );
        //REPLACED BY NEW CODE BY ALBERT:
	rho1 = asin( rho1 / dsun );
	//A similar correction was applied in build_suba.c
	
        /*  use this when rho1 is in radians */
        /* rho1 = 0.0175295299447; */

	/* solar pole vector */	// Un-commented by Albert
	spol1[0] = cos(DELTApo)*cos(ALPHApo);
	spol1[1] = cos(DELTApo)*sin(ALPHApo);
	spol1[2] = sin(DELTApo);

	/* NOT IN ORIGINAL CODE */
	/*                         // Commented by Albert
	spol1[0] = 0.0;
	spol1[1] = 0.0;
	spol1[2] = 1.0;
	*/
	
	/* Calculate R12 matrix:  R12 = Rx(a3)Ry(a2)Rz(a1) */

	/* Zero y component of sun_ob */
	Rz = rotz( -atan2(sun_ob1[1], sun_ob1[0]));
	rotvmul(sob, Rz, sun_ob1);
	/* Zero z component of Rz * sunob */
	Ry = roty( -atan2(sob[2], sob[0]));
	/* Zero y component of spol */
	rotmul(&Rtmp, Ry, Rz);
	rotvmul(r3tmp, &Rtmp, spol1);   // Note this is NOT spol2 yet, it is: "Ry(b) Rz(a) spol1".
                                        // so Albert changed it to r3tmp here (and in the next line), for clarity
	                                // I did the same last Fall in build_subA.c, so I replicate it here.
	Rx = rotx( atan2(r3tmp[1], r3tmp[2]));

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
	Rz = rotz( carlong ); /* correct */
	rotmul(&R23, Rz, Ry);

	// Compute Sun_ob3:
	rotvmul(sun_ob3, &R23, sun_ob2);

	free(Rz);
	free(Ry);

	
fprintf(stderr,"polar angle: %g radians = %g deg\n",
	      pang,pang*180./((double) M_PI));
fprintf(stderr,"Carrington longitude: %g radians =  %g deg\n",
	      carlong,carlong*180./((double) M_PI));
fprintf(stderr,"spol1: [%g, %g, %g]\n",spol1[0],spol1[1],spol1[2]);
fprintf(stderr,"spol2: [%g, %g, %g]\n",spol2[0],spol2[1],spol2[2]);
      rotvmul(r3tmp, &R23, spol2);
fprintf(stderr,"spol3: [%g, %g, %g]\n",r3tmp[0],r3tmp[1],r3tmp[2]);
fprintf(stderr,"sun_ob1: [%g, %g, %g]\n",sun_ob1[0],sun_ob1[1],sun_ob1[2]);
fprintf(stderr,"sun_ob2: [%g, %g, %g]\n",sun_ob2[0],sun_ob2[1],sun_ob2[2]);
      rotvmul(r3tmp, &R23, sun_ob2);

  fprintf(stderr, "      Computed sun_ob3:  [%3.10g, %3.10g, %3.10g]\n\n",sun_ob3[0], sun_ob3[1], sun_ob3[2]);
       
	hasdata = 1;
#include "buildrow.c"
#ifdef GEOMTEST_OUT
	fclose(fid_geomtest);
#endif
}

