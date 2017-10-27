
/* this thing reads in an mk3 rectangular 512 x 512 data file 
 *  and puts the data in the same index order as the idl file
 *  be careful because the data are in the form of 2 byte ints
 */


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

void 
main( int arc, char **argv)
{

  short im1[512][512];
  float pBval[512][512];
  float image_x[512], image_y[512], rho[512][512], eta[512][512];
  int status, i , jj;
  FILE *fid;
  
  for (i = 0; i < imsize; i++){
    x_image[i] =  MK3_PIXSIZE*((float) i - MK3_CENTER_X);
    y_image[i] =  MK3_PIXSIZE*((float) i - MK3_CENTER_Y);
  }
  for (i = 0; i < imsize; i++){
    for (jj = 0; jj < imsize; jj++){
      rho[i][jj] = (float) sqrt( (double) (x_image[i]*x_image[i] + y_image[jj]*y_image[jj])  );
      eta[i][jj] = (float) atan2( (double) (- x_image[i]), (double) y_image[jj] );
    }
  }


  fid = fopen("97d031.ch0.rpB.avg","rb");

  status = fread(im1 , sizeof (short) , 512*512 , fid );

  for (i = 0 ; i < 512; i++)
    for ( jj = 0; jj < 512; jj++)
      pBval[i][jj] = (float) im1[jj][i];






  i = 161; jj = 256;

  fprintf(stdout,"%g\n",pBval[i][jj]);

}
