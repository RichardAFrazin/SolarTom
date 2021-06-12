/* the various build definitions need to be set in builA_params.h */
 
/*====  DEFINES FOR SOLVE AND FRIENDS ====*/

#define NMATS 2                         /* total number of matrices: A1, possibly A2, plus Reg1,... */  
#define NUMBER_OF_DATA_MATRICES 1	/* number observation matrices */
#define START_TOL 1.e-5	                /* beginning iteration tolerance */
#define CHANGETOL_FACTOR 5.0	        /* divide tolerance by this factor once it's been reached */
#define FRACTIONAL_CHANGE_TOL 0.01	/* mean fractional object differenrce required to exit */
#define MINIMUM_VALUE 0 	        /* min value assigned to a pixel intensity */
#define ITMAX 5000		        /* max no. calls to minimizer fcn */
#define NSUBIT 1000		        /* number of suberiteration per call only used for CG algorithm */

#define PRINT_FILE_INFO		        /* prints file info */

// If running with NMATS=2 and NUMBER_OF_DATA_MATRICES=1,
// LAMBDA, FILESTRO (data), FILESTR1 (reg matrix) can be specified in the calling sequence.
// In that case any values provided below are overwritten.
//

// If running with NMATS > 2 then all inputs must be specified below and the calling sequence has no parameters.

#if (defined EITBUILD || defined EUVIBUILD || defined AIABUILD || defined WISPRIBUILD || defined WISPROBUILD || defined KCORBUILD || defined COMPBUILD || defined METISVLBUILD)
#define LAMBDA  {1.,1.e-6}             // LAMBDA and HUBER_FLAG should have NMATS elements. Extra elements are ignored.
#define HUBER_FLAG {0,0}
#define FILESTR0 "xxx"
#define FILESTR1 "r3_18_60_120"
//#define FILESTR0 "wisprI.Synth.CR2082.UnifLong.SciOrb12.bf4"
//#define FILESTR1 "wisprO.Synth.CR2082.UnifLong.SciOrb12.bf4"
//#define FILESTR2 "d2r_50_90_180"      // Must always be specified.
//#define FILESTR0 "wisprI.512.CircularOrbit01.60images"                        // A_outfile of first A matrix
//#define FILESTR1 "wisprO.512.CircularOrbit01.60images"                        // A_outfile of second A matrix, or first Reg matrix
//#define FILESTR2 "hlaplac_100_90_180"      // Must always be specified.
//#define FILESTR2 "identity_100_90_180"      // Must always be specified.
//#define FILESTR2 "d2r_100_90_180"      // Must always be specified.
//#define FILESTR3 "d2theta_100_90_180"      // Must always be specified.
//#define FILESTR4 "d2phi_100_90_180"      // Must always be specified.
#define MAIN_X_INFILE  "x_AWSOM_CR2082_sphere_WISPR.dat"
#define MAIN_X_OUTFILE "x_wisprIO.512.CR2082.UnifLong.SciOrb12.bf4_r3_l1e-6"

#elif (defined C2BUILD || defined CORBUILD || defined C3BUILD)
#define HUBER_FLAG {0, 0}
#define FILESTR0 "xxx"
#define FILESTR1 "r3_60_60_120"
#define MAIN_X_INFILE  "x_AWSOM_CR2081run5_WISPR_sphere_2.dat"
#define MAIN_X_OUTFILE " "
#define LAMBDA  { 1.0 , 1.e-6}

#else
#error No build that is currently understood specified
#endif

/* Cross Validation stuff (see auto_cv(_brent), cv_brent_fixed, cv_calc */
/*#define PRESERVE_SOL */   /* define to keep the solution computed at each iteration (cvcalc) */
#define CV_X_OUTFILE    "x_comp1079.dynamics.Dt2.CR2198.bf2.ri1.00.ro1.50_50_90_180_hlaplac-d2r_r1w1" /* cv solution (code adds _auto_cv suffix) - not used with -o opiton in auto_cv */
#define CV_X_INFILE     "x_comp1079.dynamics.Dt2.CR2198.bf2.ri1.00.ro1.50_50_90_180_hlaplac-d2r_r1w1" /* cv initial solution - not used with -i opiton in auto_cv */
#define AMOEBA_ITMAX 12			/* max number of fcn evals */
#define MIN_LAMBDA 1.e-7    /* min and max lambda (used by brent) */
#define MAX_LAMBDA 1.e-3  
#define BRENT_TOL  0.01*MIN_LAMBDA /* brent tollerance (in x) for stopping */
#define AMOEBA_LAMBDA {{7.e-5},{1.e-4}}   /* only used for amoeba */

/*#define AMOEBA_LAMBDA { {1.05e-7, 2.46e-6, 2.22e-7},	\
                        {4.4e-7, 1.4e-5, 8.5e-6},\
                        {1.4e-7, 2.6e-6, 2.3e-7},\
                        {3.e-9, 9.4e-6, 5.e-6} } */
