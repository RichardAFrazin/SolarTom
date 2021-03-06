/* the various build definitions need to be set in builA_params.h */



/*====  DEFINES FOR SOLVE AND FRIENDS ====*/

#define NMATS 2                      /* total number of matrices */
#define NUMBER_OF_DATA_MATRICES 1	/* number observation matrices */
#define START_TOL 1.e-5	                /* beginning iteration tolerance */
#define CHANGETOL_FACTOR 5.0	        /* divide tolerance by this factor once it's been reached */
#define FRACTIONAL_CHANGE_TOL 0.01	/* mean fractional object differenrce required to exit */
#define MINIMUM_VALUE 0 	/* min value assigned to a pixel intensity */
#define ITMAX 5000		/* max no. calls to minimizer fcn */
#define NSUBIT 600		/* number of suberiteration per call only used for CG algorithm */

#define PRINT_FILE_INFO		/* prints file info */


#if (defined EITBUILD || defined EUVIBUILD || defined AIABUILD)
#define LAMBDA  { 1.0 , 1.0, 100.0}
#define HUBER_FLAG {0, 0, 0}
#define FILESTR0 "euviA.171.fullcr2069" 
#define FILESTR2 "" 
#define FILESTR1 "hlaplac_20_90"
#define MAIN_X_INFILE  " "
#define MAIN_X_OUTFILE " "

#elif (defined C2BUILD || defined CORBUILD || defined C3BUILD)
#define HUBER_FLAG {0, 0}
#define FILESTR0 A_OUTFILE " "
#define FILESTR1 "hlaplac_60_60"
#define MAIN_X_INFILE   "x_zero_100_180"
#define MAIN_X_OUTFILE " "
#define LAMBDA  { 1.0 , 1.e-6}

#else
#error No build that is currently understood specified
#endif


/* Cross Validation stuff (see auto_cv(_brent), cv_brent_fixed, cv_calc */
/*#define PRESERVE_SOL */   /* define to keep the solution computed at each iteration (cvcalc) */
#define CV_X_OUTFILE    "x_euviA_171_200701" /* cv solution (code adds _auto_cv suffix) - not used with -o opiton in auto_cv */
#define CV_X_INFILE     "x_euviA_171_200701l.01" /* cv initial solution - not used with -i opiton in auto_cv */
#define AMOEBA_ITMAX 12			/* max number of fcn evals */
#define MIN_LAMBDA 1.e-6  /* min and max lambda (used by brent) */
#define MAX_LAMBDA 1.e-4
#define BRENT_TOL (0.01*MAX_LAMBDA) /* brent tollerance (in x) for stopping */
#define AMOEBA_LAMBDA {{7.e-5},{1.e-4}}   /* only used for amoeba */

/*#define AMOEBA_LAMBDA { {1.05e-7, 2.46e-6, 2.22e-7},	\
                        {4.4e-7, 1.4e-5, 8.5e-6},\
                        {1.4e-7, 2.6e-6, 2.3e-7},\
                        {3.e-9, 9.4e-6, 5.e-6} } */
