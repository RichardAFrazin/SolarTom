function [x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename,path)
% function [x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename, [path])

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;
 
if (type ~= 'CYLINDRICAL')
  error(['Recontruction must be on a cylindircal grid to use this ' ...
         'function']);
end

nc3 = nrad*nphi*nz;

if (nargin == 1)
  datadir = ['../', bindir];
elseif (nargin == 2)
  datadir = [path, '/'];
else
  error('assertion failed');
end

fname = strcat(datadir,filename);
disp(['filename = ',fname]);
fid = fopen(fname,'rb');

[x, cnt] = fread(fid, nc3,'float32');

fclose(fid);

if (cnt ~= nc3)
  disp(['error in fread! cnt = ',num2str(cnt),' nbins = ',num2str(nc3)]);
  stop
end

x = reshape(x,nrad,nphi,nz);

% Always compare to ind2cyl.m for correct definition

r = (grid_rmax / nrad / 2) * (2 * (1:nrad) - 1);
phi = 180 / pi * (2 * pi / nphi / 2) * (2 * (1:nphi) - 1);
z = (2 * grid_rmax / nz / 2) * (2 * (1:nz) - 1) - grid_rmax;