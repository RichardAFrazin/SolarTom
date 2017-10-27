function gen_constraint(path)
% function gen_constraint([path])

% Maximum fast wind proton flux at 1 AU.  Units flux / cm^3.  Found
% from in-situ spacecraft measurement.

max_fw_pf = 3;


% Assuming the solar wind to have an upper limit of 800 km/s and
% constant.

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

% Always compare to ind2cyl.m for correct definition

r = (grid_rmax / nrad / 2) * (2 * (1:nrad) - 1);
z = (2 * grid_rmax / nz / 2) * (2 * (1:nz) - 1) - grid_rmax;

x = zeros(nrad, nphi, nz);

constraint = zeros(nrad, nz);

for i=1:nrad
  for j=1:nz
    r_sphere = sqrt(r(i)^2 + z(j)^2);
    constraint(i,j) = 3*(215/r_sphere)^2;  % lower bound on electron density
  end
end

for i=1:nphi
  x(:,i,:) = constraint;
end

fname = ['cstr', num2str(nrad), num2str(nphi), num2str(nz)];

if (nargin == 1)
  fname = [path, fname];
else
  fname = ['../', bindir, fname];
end

disp(['Writing constraint file: ', fname]);

fid = my_fopen(fname, 'w');

count = fwrite(fid, x(:), 'float32');

if (count ~= nrad*nphi*nz)
  error(['Error writing to file: ', fname]);
end

fclose(fid);

disp('done');