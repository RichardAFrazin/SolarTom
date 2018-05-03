function gen_bdy(path)
% function gen_bdy([path])

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

fname = ['bdy', num2str(nrad), num2str(nphi), num2str(nz), '.' ...
         num2str(round(rmin*100)), '.', num2str(round(rmax*100))];

if (nargin == 1)
  fname = [path, fname];
else
  fname = ['../', bindir, fname];
end

fid = fopen(fname, 'w');
if (fid == -1)
  error(['Failed to open file for writing: ', fname]);
end

r = linspace(grid_rmax/nrad,grid_rmax,nrad);
phi = linspace(0,360*(nphi-1)/nphi,nphi); 
z = linspace(-grid_rmax + (2*grid_rmax)/nz,grid_rmax,nz);

nc3 = nrad * nphi * nz;

b = zeros(nrad,nphi,nz);

mask = zeros(nrad, nz);


for r_i=1:nrad
  for z_i=1:nz
    dist = sqrt(z(z_i)^2 + r(r_i)^2);
    if (dist > rmax )%| dist < rmin - 3*rmax/nrad )
      mask(r_i,z_i) = 0;
    else
      mask(r_i,z_i) = 1;
    end
  end
end

for i=1:nphi
  b(:,i,:) = mask;
end

b = b(:);

count = fwrite(fid, b, 'int32');

if (count ~= nc3)
  error(['Error writing to file: ', fname]);
end

fclose(fid);