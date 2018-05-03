function [i_rad,i_phi,i_z] = cyl2ind(rad, phi, z)
% function [i_rad,i_phi,i_z] = cyl2ind(rad, phi, z))
%
% Convert from heliocentric cylindrical coordinates to data
% cylinder indices.

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

deltagrid = 2 * grid_rmax / nz;

if (z == grid_rmax)
  i_z = nz;
else
  i_z = floor((z + grid_rmax) / deltagrid) + 1;
end

if (rad == grid_rmax)
  i_rad = nrad;
else
  i_rad = floor((rad / grid_rmax) * nrad) + 1;
end

if (phi == 2 * pi)
  i_phi = n_phi;
else
  i_phi = floor((phi / 2 / pi) * nphi) + 1;
end