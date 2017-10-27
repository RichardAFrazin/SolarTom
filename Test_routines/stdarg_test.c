/*   stdarg_test.h by Richard Frazin 5/2008
 * This tests my understanding of functions with variable
 * input arguments.
 *
 *  It is necessary to inform the function, in some way,
 *    of how many variables to expect.
 *
 * to complile: gcc -o stdarg_test stdarg_test.c -lm 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdarg.h>


/* this is a function from the Gnu C programming guide 
 *   count contains the number of arguments */
/*
int add_em_up (int count,...)
{
  va_list ap;
  int i, sum;

  va_start (ap, count);     
  sum = 0;
  for (i = 0; i < count; i++)
    sum += va_arg (ap, int);   

  va_end (ap);               
  return sum;
}

int
main (void)
{
  printf ("%d\n", add_em_up (3, 5, 5, 6));
  printf ("%d\n", add_em_up (10, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10));
  return 0;
}
*/

int stdarg_test(char *, int, int, ...);

int main(int argc, char **argv){
  char foo1[128], foo2[128], foo3[128];
  int  bar, result;

  if ((argc < 3) || (argc > 5)){
    fprintf(stdout,"Use 2,3 or 4 command line arguments: %s <string1> <integer> [<string2>] [<string3>] \n",argv[0]);
    return(1);
  }

  fprintf(stdout,"argc = %d\n",argc); fflush(stdout);
  if (argc == 3){
    strcpy(foo1,argv[1]);
    bar = atoi(argv[2]);
    result = stdarg_test(foo1, bar, 0);
    fprintf(stdout,"ans = %d, %s\n",result,foo1);
    return(0);
  } else if (argc == 4){
    strcpy(foo1,argv[1]);
    strcpy(foo2,argv[3]);
    bar = atoi(argv[2]);
    result = stdarg_test(foo1, bar, 1, foo2);
    fprintf(stdout,"ans = %d, %s\n",result,foo1);
    return(0);
  } else if (argc == 5){
    strcpy(foo1,argv[1]);
    strcpy(foo2,argv[3]);
    strcpy(foo3,argv[4]);
    bar = atoi(argv[2]);
    result = stdarg_test(foo1, bar, 2, foo2, foo3);
    fprintf(stdout,"ans = %d, %s\n",result,foo1);
    return(0);
  }  

}

/* noptarg is the number of optional arguments */
int stdarg_test(char *cc, int jj, int noptarg, ...){
  char *dd;
  int kk, count;
  va_list ap;

  strcat(cc,"_X");
  kk = 5 + jj;

  if (noptarg == 0){
    return(kk);
  }

  va_start(ap,noptarg);
  for (count = 0; count < noptarg; count++){
    dd = va_arg(ap, char *);
    strcat(cc,"_");
    strcat(cc,dd);
  }
  va_end(ap);

  return(kk);
}

