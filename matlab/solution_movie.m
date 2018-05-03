function [x,rad,lat,lon] = solution_movie(filename,path,nrad,ntheta,Rmin,Rmax,endian);
%function [x,rad,lat,lon] = solution_movie(filename,path,nrad,ntheta,Rmin,Rmax,endian);
%
% This calls readtom_sph.m and makes a movie showing the solution at 
%    constant radius.
%
%x = output array (nrad,ntheta,2*ntheta)
%r = radial grid  linspace(1+(Rmax-1)/nrad,Rmax,nrad)
%lat = latitude grid in deg (ntheta)
%lon = longitude grid in deg (2*ntheta)
%filename - duh?
%path - include final slash
%nrad = # of radial grid points
%ntheta = # of phi grid points
%Rmax = max radius of computation grid in Rsun
%endian = an optional argument for endian 
%    or other optional format string. usually 'ieee-be' or 'ieee-le'
%
%

if (nargin == 7)
  endianstring = endian;
else 
   endianstring = 'ieee-le';
end

[x,rad,lat,lon] = readtom_sph(filename,path,nrad,ntheta,Rmin,Rmax,endianstring);

if (min(min(min(x))) < 0)
   disp('Warning: non-positive solution.  linear scale displayed.');
end


figure;
for k = 1:nrad
  im = reshape(x(k,:,:),ntheta,2*ntheta);
  
  if (min(min(min(x))) < 0)
    imagesc(lon,lat,im); colorbar;
    title(['Ne  r(',num2str(k),') = ',num2str(rad(k))]);
  else
    imagesc(lon,lat,sqrt(im)); colorbar;
    title(['SQRT(Ne)  r(',num2str(k),') = ',num2str(rad(k))]);
  end
  set(gca,'YDir','normal','FontSize',14);
  pause;
end

return;
