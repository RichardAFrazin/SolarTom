function [x,rad,lat,lon,time] = dynamic_readtom_sph(filename,baseproj,path,nrad,ntheta,Rmin,Rmax,endian)
%function [x,rad,lat,lon,time] = dynamic_readtom_sph(filename,path,nrad,ntheta,Rmin,Rmax,endian)
%
%This reads the tomographic output for the hollow sphere geometry
%
%x = output array (nrad,ntheta,2*ntheta)
%r = radial grid  linspace(1+(Rmax-1)/nrad,Rmax,nrad)
%lat = latitude grid in deg (ntheta)
%lon = longitude grid in deg (2*ntheta)
%time = time in days starting at 0
%filename - dynamic solution filename
%baseproj - filename base of projection matrix
%path - include final slash
%nrad = # of radial grid points
%ntheta = # of phi grid points
%Rmax = max radius of computation grid in Rsun
%endian = an optional argument for endian 
%    or other optional format string. usually 'ieee-be' or 'ieee-le'

if (nargin == 8)
  endianstring = endian;
else 
   endianstring = 'ieee-le';
end

nphi = 2*ntheta;
nc3 = nrad*ntheta*nphi;

fname = strcat(path,filename);
disp(['filename = ',fname]);
fid = fopen(fname,'rb');
if (fid < 0 ) 
    disp('dynamic_readtom_sph: bad file name');
    return
end

[x, cnt] = fread(fid, inf,'float32',endianstring);
fclose(fid);

ntime = cnt/nc3;
if ( (fix(ntime) ~= ntime) | (cnt <= 0) ) 
    disp('dynamic_readtom_sph: Problem reading x, bad file length');
    return;
end

x = reshape(x,nrad,ntheta,nphi,ntime);

nphi = 2*ntheta;
% these are the midpoints of the cells
rad = linspace(Rmin + (Rmax-Rmin)/nrad/2, Rmax - (Rmax-Rmin)/nrad/2,nrad);
lat = linspace(-90 + 90/ntheta,90 - 90/ntheta,ntheta);
lon = linspace(360/nphi/2, 360 - 360/nphi/2, nphi );

%read date file
mjd = binfileread([path,'date_'],baseproj,'float64');
if ((mjd(1) == -1) )
    disp('dynamic_readtom_sph: problem reading:')
    disp([directory,'date_',baseproj]);
    return;
end
if (length(mjd) ~= ntime)
    disp('dynamic_readtom_sph: Problem reading date file, bad file length');
    disp([directory,'date_',baseproj]);
    return;
end
time=mjd-mjd(1);







