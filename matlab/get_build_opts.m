% [nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
% imsize, build, type, outfile ] = get_build_opts
%
% nrad:      (scalar) number of radial bins
% nphi:      (scalar) number of theta bins
% nz:        (scalar) number of z bins
% rmin:      (scalar) minimum radius
% rmax:      (scalar) maximum radius
% grid_rmax: (scalar) maximum radius on the computation grid
% bindir:    (string) path to bindata
% x_infile:  (string) initial x
% cv_infile: (string) initial x for cv
% imsize:    (scalar) input image is imsize x imsize
% build:     (string) the build string
% type:      (string) build type, CYLINDRICAL or CARTISIAN
% outfile:   (string) file name of the output of solve
