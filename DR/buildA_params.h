/*#define AIABUILD*/
#define EUVIBUILD
/* #define CORBUILD */
/*#define EITBUILD*/
/* #define C3BUILD */
/* #define C2BUILD */
#if (defined C2BUILD || defined C3BUILD)
/*#define NRL */      /*use for NRL        calibration */
#define MARSEILLES   /*use for Marseilles calibration */
#endif

/* differential rotation transition scale-height */
#define DIFFROTPARAMVALUE 0.1

#ifdef C2BUILD
#define RMAX 8.2  /*6.1*/           /* outer radius of computation ball */
#define RMIN 2.2            /* innner radius (hollow  sphere)   */
#define NZ     130 
#define NCELLS  90	/* cartesian: object has NCELLS^3 elements */
#define NRAD   60  /*40*/
#define NTHETA 60           /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE    512	/* size of C2 images (pixels) */
#define BINFAC    3	     /* binning factor for C2 images (pixels) */
#define DELTA     0.0	/* delta vector */
#define INSTR_RMIN      2.2
#define INSTR_RMAX      6.1
#define PIXSIZE  23.8	/* arcsec per pixel */
#ifdef NRL
typedef float PB_IMTYPE;
#endif
#ifdef MARSEILLES
typedef double PB_IMTYPE;
#endif
#define DATADIR    TOMROOT"DATA/c2Mars/1996/"
#define CONFSTRING DATADIR"list.c2mars.pb.1996.0801.0814.txt"
#define A_OUTFILE     "c2mars.pb.1996.0801.0814-2.2-8.2nr60nt60" /* suffix of A matrix ouput files */

#elif defined C3BUILD
#define RMAX 10.2            /* outer radius of computation ball */
#define RMIN 4.2            /* innner radius (hollow  sphere)   */
#define NZ     130 
#define NCELLS  90	/* cartesian: object has NCELLS^3 elements */
#define NRAD   60  
#define NTHETA 60           /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE    512	/* size of C2 images (pixels) */
#define BINFAC    1	     /* binning factor for C2 images (pixels) */
#define DELTA     0.0	/* delta vector */
#define INSTR_RMIN      4.1
#define INSTR_RMAX      18.0
#define PIXSIZE  (27.44*1024/IMSIZE)  /* arcsec per pixel */
#ifdef NRL
typedef float PB_IMTYPE;
#endif
#ifdef MARSEILLES
typedef double PB_IMTYPE;
#endif
#define DATADIR    TOMROOT"DATA/LASCO/2006/"
#define CONFSTRING DATADIR"testc3_060806.txt"
#define A_OUTFILE     "testc3_060806.txt" /* suffix of A matrix ouput files */

#elif defined CORBUILD  /* for COR1 and COR2(?)*/
#define NRL /* pB scaling */
#define RMAX 6.1            /* outer radius of computation ball */
#define RMIN 1.6            /* innner radius (hollow  sphere)   */
#define NZ     130 
#define NCELLS  90	/* cartesian: object has NCELLS^3 elements */
#define NRAD   60  
#define NTHETA 60           /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE    512	/* size of COR1 images (pixels) */
#define BINFAC    2	     /* binning factor for images (pixels) */
#define DELTA     0.0	/* delta vector */
#define INSTR_RMIN      1.6
#define INSTR_RMAX      4.2
#define PIXSIZE  (7.5*1024.0/IMSIZE)	/* COR1 arcsec per pixel */
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/COR1/cr2106/"
#define CONFSTRING DATADIR"test_B.txt"
#define A_OUTFILE     "test_B" /* suffix of A matrix ouput files */



/*   FOR EUVI, EIT Only! */
/* With this scaling (v,w in units of RSUN, y units of modified DN/s),
 *  after x is determined it needs to be DIVIDED by 
 *  RSUN*1.e-12 = 0.0696   */

#elif defined EITBUILD
#define RMAX 1.26           /* outer radius of computation ball */
#define RMIN 1.0            /* innner radius (hollow  sphere)   */ 
#define NRAD 13             /* radial bins */
#define NTHETA 60           /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE 512      /* this operates on processed files rebinned to 512 */
#define DELTA 0.0
#define INSTR_RMAX 1.26
#define PIXSIZE 5.26    /* arcsec per pixel (512x512) */
#define BINFAC 2
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/EIT/2008/"
#define CONFSTRING DATADIR"list.txt"
#define A_OUTFILE      "crap"

#elif defined EUVIBUILD
/* don't accept data between INNER_REJECT_RAD and OUTER_REJECT_RAD (due to optical depth) */
#define RING_REJECT
#define RMAX 1.205
#define RMIN 1.005
#ifdef  RING_REJECT
#define INNER_REJECT_RAD 0.98   /* 0.98 */ 
#define OUTER_REJECT_RAD 1.025 /* 1.025 */
#endif
#define NRAD 20
#define NTHETA 90
#define NPHI (NTHETA * 2)
#define IMSIZE 1024
#define DELTA 0.0
#define INSTR_RMAX RMAX
#define PIXSIZE (1.589*2048./IMSIZE) /*A is 1.588 ''/pix, B is 1.590 in 2048 mode*/
#define BINFAC 4
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/EUVI/CR2084/A171/" 
#define CONFSTRING DATADIR"list.A171.b4.little.txt"        
#define A_OUTFILE  "CR2084.284A.20.90_DR0.1e"


#elif defined AIABUILD
#define RING_REJECT  /* don't accept data between INNER_REJECT_RAD and OUTER_REJECT_RAD (due to optical depth) */
#define RMAX 1.26
#define RMIN 1.0
#ifdef RING_REJECT 
#define INNER_REJECT_RAD 0.
#define OUTER_REJECT_RAD 1.015
#endif
#define NRAD 26
#define NTHETA 90 
#define NPHI (NTHETA * 2)
#define IMSIZE 1024
#define DELTA 0.0
#define INSTR_RMAX 1.26
#define PIXSIZE (0.6*4096./IMSIZE)  
#define BINFAC 4
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/AIA/094/" 
#define CONFSTRING DATADIR"list.094.processed.binned.txt"        
#define A_OUTFILE  "AIA.CR2106.094.nr26.irm1.26"

#else
#error No build that is currently understood specified
#endif
