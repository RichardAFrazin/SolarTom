function [x,rad,lat,lon] = readtom_sph(filename,path,nrad,ntheta,Rmin,Rmax,endian)
%function [x,rad,lat,lon] = readtom_sph(filename,path,nrad,ntheta,Rmin,Rmax,endian)
%
%This reads the tomographic output for the hollow sphere geometry
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

if (nargin == 7)
  endianstring = endian;
else 
   endianstring = 'ieee-le';
end

nphi = 2*ntheta;
nc3 = nrad*ntheta*nphi;

fname = strcat(path,filename);
disp(['filename = ',fname]);
fid = fopen(fname,'rb');

[x, cnt] = fread(fid, nc3,'float32',endianstring);

fclose(fid);

if (cnt ~= nc3)
  disp(['error in fread! cnt = ',num2str(cnt),' nbins = ',num2str(nc3)]);
  return;
end

x = reshape(x,nrad,ntheta,nphi);

nphi = 2*ntheta;
% these are the midpoints of the cells
rad = linspace(Rmin + (Rmax-Rmin)/nrad/2, Rmax - (Rmax-Rmin)/nrad/2,nrad);
lat = linspace(-90 + 90/ntheta,90 - 90/ntheta,ntheta);
lon = linspace(360/nphi/2, 360 - 360/nphi/2, nphi );



