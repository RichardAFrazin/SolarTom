function [rad,phi,z] = ind2cyl(i_rad, i_phi, i_z)
% function [rad,phi,z] = ind2cyl(i_rad, i_phi, i_z)
%
% Convert from data cylinder indices to heliocentric cylindrical
% coordinates.

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

z = (2 * grid_rmax / nz / 2) * (2 * i_z - 1) - grid_rmax;

rad = (grid_rmax / nrad / 2) * (2 * i_rad - 1);

phi = (2 * pi / nphi / 2) * (2 * i_phi - 1);