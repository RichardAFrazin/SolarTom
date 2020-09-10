// #define AIABUILD
// #define EUVIBUILD
// #define CORBUILD
// #define EITBUILD
// #define C3BUILD
   #define C2BUILD 
// #define WISPRIBUILD
// #define WISPROBUILD
// #define KCORBUILD
// #define COMPBUILD

/* Not using this for now
   #if (defined WISPIRIBUILD || defined WISPROBUILD)
   #define Orb_1  // Select Orbit Number, add suffix with right orbit number after "_"
   #endif
*/

//#define TESTBAND // Set to make compare.c compute only a full "band"
                 // of horizontal pixels of BAND_WIDTH_PX pixels around CRPIX2.
#if defined TESTBAND
#define BAND_WIDTH_PX 100
#endif

// #define NONUNIFORMRAD // Set if radial grid is not uniform.

#ifdef NONUNIFORMRAD
#define GRID_FILENAME "non_uniform_grid.txt"
#else
#define GRID_FILENAME "uniform_grid.txt"
#endif
   
#if (defined C2BUILD || defined C3BUILD)
// #define NRL            // use for NRL        calibration 
#define MARSEILLES     // use for Marseilles calibration
#endif

#ifdef WISPRIBUILD
#define RMIN  2.0             /* RMIN reached by the grid */
#define RMAX  214.5        //14.//14.196949// 214.5            /* RMAX reached by the grid */
#define NZ     130
#define NCELLS  90	     /* cartesian: object has NCELLS^3 elements */
#define NRAD    100         //48//20//100
#define NTHETA  90           /* polar angle bins */
#define NPHI (NTHETA * 2)    /* azimuthal angle bins */
#define IMSIZE    512	     /* size of WISPR images (pixels), expanded 1920->2048 in height to make them square in first version*/
#define BINFAC    4	     /* binning factor for WISPRI images (pixels) */
#define DELTA     0.0	     /* delta vector */
#define INSTR_RMIN  2.0    /* Set as the range of radii over the whole image series over a full 0.5AU->0.5AU orbit*/
#define INSTR_RMAX  1000.     //1000. // 90.    
#define PIXSIZE     (66.621094*2048/IMSIZE)  /* arcsec per pixel */
typedef float PB_IMTYPE;
#define DATADIR     TOMROOT"DATA/wisprI/CR2082_UnifLong/"
#define CONFSTRING  DATADIR"list.wisprI.Synth.CR2082.UnifLong.SciOrb12.txt"
#define A_OUTFILE               "wisprI.Synth.CR2082.UnifLong.SciOrb12.bf4"     /* suffix of A matrix ouput files */

#elif defined WISPROBUILD
#define RMIN 2.0             /* RMIN and RMAX set as in WISPRIBUILD, see notes above */
#define RMAX 214.5     //14.//14.196949//214.5            
#define NZ     130 
#define NCELLS  90   	     /* cartesian: object has NCELLS^3 elements */
#define NRAD    100            //48//20//100 
#define NTHETA  90           /* polar angle bins */
#define NPHI (NTHETA * 2)    /* azimuthal angle bins */
#define IMSIZE    512	     /* size of WISPR images (pixels), expanded 1920->2048 in height to make them square in first version */
#define BINFAC    4	     /* binning factor for WISPRO images (pixels) */
#define DELTA     0.0	     /* delta vector */
#define INSTR_RMIN   2.0       //7.0
#define INSTR_RMAX   1000.   //1000.   //110. 
#define PIXSIZE     (96.679680*2048/IMSIZE)  /* arcsec per pixel */
typedef float PB_IMTYPE;
#define DATADIR     TOMROOT"DATA/wisprO/CR2082_UnifLong/"
#define CONFSTRING  DATADIR"list.wisprO.Synth.CR2082.UnifLong.SciOrb12.txt"
#define A_OUTFILE               "wisprO.Synth.CR2082.UnifLong.SciOrb12.bf4"     /* suffix of A matrix ouput files */

#elif defined KCORBUILD
#define RMIN   1.05          /* innner radius (hollow  sphere)   */
#define RMAX   2.25          /* outer radius of computation ball */
#define NRAD    120
#define NTHETA   90         /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE  1024	    /* size of C2 images (pixels) */
#define BINFAC    2	    /* binning factor for C2 images (pixels) */
#define DELTA     0.0	    /* delta vector */
#define INSTR_RMIN      1.09
#define INSTR_RMAX      2.00
#define PIXSIZE     (5.643*1024/IMSIZE)     /* arcsec per pixel */
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/kcor/CR2198/10MinAvg/"
#define CONFSTRING DATADIR"list_new_reduced.txt"
#define A_OUTFILE         "KCOR.CR2198.13imgs-reduced.bf2.ri1.05.ro2.25_Inst_1.09_2.00_120_90_180_dropneg" /* suffix of A matrix ouput files */

#elif defined COMPBUILD
#define RMIN   1.00         /* innner radius (hollow  sphere)   */
#define RMAX   1.50         /* outer radius of computation ball */
#define NRAD     50
#define NTHETA   90         /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define BINFAC    2	    /* binning factor for C2 images (pixels) */
#define DELTA     0.0	    /* delta vector */
#define INSTR_RMIN      1.09
#define INSTR_RMAX      1.30
#define IMSIZE       620	    /* size of COMP images (pixels) */
#define PIXSIZE     (4.350*620/IMSIZE)     /* arcsec per pixel */
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/comp/1074/CR2198/"
#define CONFSTRING DATADIR"list_total_intensity_Dt2.txt"
#define A_OUTFILE         "comp1074.dynamics.Dt2_CR2198.bf2.ri1.00.ro1.50_50_90_180_INSTR_RMIN_1.09_ABBEY" /* suffix of A matrix ouput files */

#elif defined C2BUILD
#define RMIN        2.0            /* innner radius (hollow  sphere)   */
#define RMAX        8.0            /* outer radius of computation ball */
#define INSTR_RMIN  2.5
#define INSTR_RMAX  6.0
//#define NZ     130
//#define NCELLS  90	/* cartesian: object has NCELLS^3 elements */
#define NRAD        60  
#define NTHETA      90      /* polar angle bins */
#define NPHI (NTHETA * 2)   /* azimuthal angle bins */
#define IMSIZE    512	    /* size of C2 images (pixels) */
#define BINFAC    4	    /* binning factor for C2 images (pixels) */
#define DELTA     0.0	    /* delta vector */
#define PIXSIZE  (23.8*512/IMSIZE)   /* arcsec per pixel */
#ifdef NRL
typedef float PB_IMTYPE;
#endif
#ifdef MARSEILLES
typedef double PB_IMTYPE;
#endif
#define DATADIR    TOMROOT"DATA/c2/CR2219/"
#define CONFSTRING DATADIR"list_prep.txt"
#define A_OUTFILE         "LASCO-C2-pB_56imgs_Rmin2.0_Rmax8.0_IRmin2.5_IRmax6.0_60x90x180_BF4" /* suffix of A matrix ouput files */

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
#define DATADIR    TOMROOT"DATA/C2MARS/2008/"
#define CONFSTRING DATADIR"testc3_060806.txt"
#define A_OUTFILE         "testc3_060806.txt" /* suffix of A matrix ouput files */

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
#define A_OUTFILE         "test_B" /* suffix of A matrix ouput files */

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
#define A_OUTFILE         "crap"

#elif defined EUVIBUILD
/* don't accept data between INNER_REJECT_RAD and OUTER_REJECT_RAD (due to optical depth) */
#define RING_REJECT
#define INSTR_RMAX 1.5
#define RMAX       2.0
#define RMIN       1.0
#ifdef  RING_REJECT
#define INNER_REJECT_RAD 0.98  /* 0.98 */ 
#define OUTER_REJECT_RAD 1.02 /* 1.025 */
#endif
#define NRAD    100
#define NTHETA  90
#define NPHI (NTHETA * 2)
#define IMSIZE 1024
#define DELTA 0.0
#define PIXSIZE (1.589*2048./IMSIZE) /*A is 1.588 ''/pix, B is 1.590 in 2048 mode*/
#define BINFAC 4
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/euvi/CR2081/A284/"
#define CONFSTRING DATADIR"list.A284.b4.nodecon.test"
#define A_OUTFILE         "euvi_A284_b4_nodecon_100_90_180_Rmax2.0_InstrRmax1.5_bf4"
			  
#elif defined AIABUILD
#define RING_REJECT  /* don't accept data between INNER_REJECT_RAD and OUTER_REJECT_RAD (due to optical depth) */
#define NRAD         26
#define RMAX       1.26
#define RMIN       1.00
#define INSTR_RMAX 1.26
#ifdef RING_REJECT 
#define INNER_REJECT_RAD 0.00
#define OUTER_REJECT_RAD 1.02
#endif
#define NTHETA 90 
#define NPHI   (NTHETA * 2)
#define IMSIZE 1024
#define DELTA 0.0
#define PIXSIZE (0.6*4096./IMSIZE)  
#define BINFAC 4
typedef float PB_IMTYPE;
#define DATADIR    TOMROOT"DATA/aia/CR2209/171/" 
#define CONFSTRING DATADIR"list.171.processed.selected_full.b_H1"
#define A_OUTFILE         "AIA_CR2209_171b_H1_Rmin1.00_Rmax1.30_Nr30_InstRmax1.26_bf4_TEST"

#else
#error No build that is currently understood specified
#endif
