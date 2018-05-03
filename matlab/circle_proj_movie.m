function [n_images,datacube] = ciricle_proj_movie(r0,matrix_fn,x_fn)
%function [n_images,datacube] = ciricle_proj_movie(r0,matrix_fn,x_fn)
%
% warning: make sure the correct grid is used in circle_proj_sph
%
% n_images - no. of images in log file
% datacube - 3D array containing projection movie
% r0 - radius of shell (units of Rs)
% matrix_fn - matrix suffix name 
% x_fn - filename of solution

mfn = ['/Users/frazin/tomography/bindata/','log_',matrix_fn]
fid = fopen(mfn,'rt');

datacube = [];
imsize = 300;
n_images = 0;

while (feof(fid) == 0)
  tline = fgetl(fid);
  datestamp = tline(1:13);
  start = findstr('cl=',tline) + 3;
  stop =  findstr('deg, polar_ang',tline)-2;
  carlong = str2num(tline(start:stop));
  imname = tline(1:start-2); %filename from log file
  start = findstr('so1=[',tline) + 4;
  stop =findstr('] Rs',tline);
  sun_ob1 = str2num(tline(start:stop));
  
  cproj = circle_proj_sph(x_fn, r0, imsize, carlong, sun_ob1);
  fname_out = ['proj_','_',x_fn,'_',matrix_fn,'_',datestamp,'_r',num2str(r0)];
  disp(fname_out);
  fido = fopen(['/Users/frazin/tomography/bindata/Compare/',fname_out],'wb');
  fwrite(fido,cproj,'float32');
  fclose(fido);
  n_images = n_images+1;
  datacube = cat(3,datacube,cproj);
  
end
fclose(fid);
return;