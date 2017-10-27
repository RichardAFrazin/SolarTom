/*
 *=======================================================================
 *			Verify FMINBR routine
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <assert.h>


#define EPSILON 1.e-8

double fminbr(double a, double b, double (*f)(), double tol);
double fminbr_3point(double, double, double, double, double, double (*f)(), double);



static int counter;			/* Iteration counter	*/

test(a,b,f,msg)			/* Run a test			*/
double a,b;				/* Range the root is seeked for */
double (*f)(double x);			/* Functiom under examination	*/
char * msg;				/* Explanation message		*/
{
  double minloc;
  counter = 0;
  printf("\nFor function %s\nin [%g,%g] min found is at\t%.9e\n",msg,a,b,
	 (minloc=fminbr(a,b,f,EPSILON)) );
  printf("Min function value found\t%.4e\nNo. of iterations\t\t%d\n",
	 (*f)(minloc), counter);
}

double f1(x)				/* Test of the Forsythe book	*/
double x;
{
  counter++;
  return (pow(x,2)-2.0)*x - 5.0;
}

double f2(x)				/* Modified test 1            	*/
double x;
{
  double val;
 
  counter++;
  val = pow( (pow(x,2)-2.0)*x - 5.0, 2 );
  fprintf(stdout,"(%g: %g) ",x,val);
  fflush(stdout);
  return(val);
}

double f3(x)				/* My test                  	*/
double x;
{
  counter++;
  return pow( cos(x) - x,2 ) - 2;
}

double f4(x)				/* My test                  	*/
double x;
{
  counter++;
  return pow( sin(x) - x,2 ) + 1;
}



int main(int argc, char **argv){
  double (*f)(double x);
  double a,b,min,tol;
  
  a = 2.5; b = 20.; tol = 1.e-4;

  f = f2;
  /*  fminbr_3point is less effective than fminbr
  min = fminbr_3point(a,b,1.4,2.0,0.25,f,tol); 
  */
  min = fminbr(a,b,f,tol); 

  fprintf(stderr,"\nmin is at %g\n",min);

  /*min = fminbr(a,b,f,tol);
    fprintf(stderr,"min is at %g\n",min);*/

  return(0);
 }

/* original function that calls test */

not_main()
{
  test(0.0,1.0,f1,"x^3 - 2*x - 5");
  printf("Exact min is at\t\t0.81650\n");

  test(2.0,3.0,f2,"(x^3 - 2*x - 5)^2");
  printf("Exact root is \t\t2.0945514815\n");

  test(2.0,3.0,f3,"(cos(x)-x)^2 - 2");
  test(-1.0,3.0,f3,"(cos(x)-x)^2 - 2");
  test(-1.0,3.0,f4,"(sin(x)-x)^2 + 1");
}

