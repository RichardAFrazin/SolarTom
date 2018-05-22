#ifndef HEADERS_H
#define HEADERS_H
#endif

#include <stdio.h>
#include <math.h>
#include "rcs_llist.h"

//set root directory!
#include "tomroot.h"

//grid settings and data location settings
#include "buildA_params.h"

#define BINDIR  TOMROOT"bindata/"	/*location of binary input/output files */

#define DROP_NEG_PB  /* If defined, then all negative data values and their associated lines of sight are ignored */

#include "solve_cv_params.h"

#define ORBIT_URL_VERSION_MAX 9
#define ORBIT_FILE_DIR TOMROOT"Orbits/"
#define ORBIT_SCRATCH_DIR TOMROOT"Orbits/" 
#define ORBIT_URL "http://sohowww.nascom.nasa.gov/data/ancillary/orbit/predictive/"
#ifdef ORBIT_URL
#define WGETCOMMAND "/usr/bin/wget --quiet -P"  /* put the wget binary and path */
#endif


#if (defined C2BUILD || defined C3BUILD || defined CORBUILD || defined WISPRIBUILD || defined WISPROBUILD || defined KCORBUILD)
#define THOMSON		/* VdH Thomson scattering calculation */
#elif (defined EITBUILD || defined EUVIBUILD || defined AIABUILD || defined COMPBUILD)
#define RADON      /* unweighted LOS integral for EUV emissivity */
#else
#error Undefined LOS weighting!
#endif

#define NBINS (NRAD*NTHETA*NPHI)  // for hollow_sphere

/* constants */
#define ARCSECRAD (M_PI / 180. / 3600.)
#define RSUN 6.957e5		/* in km */
#define QLIMB 0.63
#define CONST 1.2497e-15	/* (3/16)*(1e10 * Thompson X-section) */
#define ALPHApo (286.13 * M_PI/180.0)	/* J2000 solar pole coords */
#define DELTApo (63.87 * M_PI/180.0)
#define MAXPATH 256		/* max string size */
#define MAX_NUMBER_OF_IMAGES 1000  



/*==== THIS IS THE MAIN DATA STRUCTURE FOR SOLVE ====*/
/*
 * iB is the vector of row indecies 
 * vB                  values
 * n[nc3+1] is the locations in the i and v arrays where 
 *    each column starts, except for the last element which
 *    has a value of 1 higher than end of that last column
 * r  is the remainder vector
 * y is the old copy of the remainder vector
 * delta is the huber function vector
 * huber is a flag: 1 robust estimate (huber fcn)
 *                  0 quadratic error penalization
 * nf = number of rows in NONSPARSE equivalent of matrix
 * file_id file identification string
 */

struct matrix {
  int *iB;
  float *vB;
  float *r, *y, *delta;
  int n[NBINS + 1];
  FILE *fid_n;
  FILE *fid_y, *fid_delta;
  FILE *fid_ind;
  FILE *fid_val;
  int huber;
  int nf;
  char file_id[MAXPATH];
};

/* sparse array stucture used for various operations
 *   ncol = # number of columns
 *   others are same as above        
 */

struct sparse {
  int nf;
  int ncol;
  int *n;
  int *iB;
  float *vB;
};


/*==== FUNCTION DECLARATIONS AND STUFF ====*/

#define MAX(x,y) ( ((x) >= (y)) ? (x) : (y) )
#define MIN(x,y) ( ((x) <= (y)) ? (x) : (y) )

typedef double Rot[3][3];

/* r3misc */
void r3norm(double *);
double r3dot(double *, double *);
void r3scalmul(double *, double);
void r3cross(double *, double *, double *);
void r3eq(double *, double *);
void r3add(double *, double *, double *);

/* build_A */

void build_subA(char *, rcs_llist *, float **, float **, int *, FILE *, FILE *);

/* grids */
double* rad_bin_boundaries(int bin);
int rad_bin_number(double dist);
void print_grid(void);
double GridDivision(double num, double denom);

void get_things_C2(char *, double *, double *);
int doublecompare(const void *, const void *);

/* rots */
Rot *rotx(double);
Rot *roty(double);
Rot *rotz(double);
void rotmul(Rot *, Rot *, Rot *);
void rotvmul(double *, Rot *, double *);

/* solve */
int solve(float *, const int *, float *, float *, float *, char *, char *, int, ...);

/* fess_hu */
void fess_hu(float *, float *, float *, struct matrix *, int *, float *);

/* cg_quad */
int cg_quad(float *, float *, int, struct matrix *);

/* normcalc */
float normcalc(char *, char *);

/* sparse */
void vmultSparse(struct sparse, float *, float *, float);
void vmult_transposeSparse(struct sparse, float *, float *, float, int);
void matmultSparse(struct sparse, struct sparse, char *);
void row_extract(struct sparse A,
                 float *y, float *delta,
                 char *fname_d, char *fname_dd,
                 int rmin, int rmax);
void row_extract_mod(struct sparse A,
                 float *y, float *delta,
		 char *fname_d, char *fname_dd,
		     int rmin, int rmax, int rstart, int rend);
void load_sparse(struct sparse *A, char *name, char *dir);
void free_sparse(struct sparse *A);
void load_vf(float **v, int n, char *fname);
void free_vf(float **v);

/* get_orbit */
void get_orbit(char *idstring, double *sun_ob, double *carlong, double *mjd);
void get_orbit_wrapper(int year, int month, int day,
                       int hour, int minute,
                       double *sun_ob, double *carlong, double *mjd);

/* amoeba */
//int amoeba(float p[NDIM + 1][NDIM], float *y, float (*func)(float *), int *iter);

/* cvcalc */
float cvcalc(double, char *, char *, char *, int, ...);

/* multiple lambda version */
/*float cvcalc(float *function_parameters,
             char *norm_matrix_name,
             char *solve_matrix_name,
             char *x_infile, char *x_outfile);*/

