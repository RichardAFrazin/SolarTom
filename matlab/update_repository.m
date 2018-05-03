function update_repository()
% function update_repository()

n_long = 128;
n_lat = 64;
r0 = 1.15;


repository = '/usr/local/arc/srt/repository';

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
  imsize, build, type, outfile ] = get_build_opts;

median_names = {
    'x_oct5_18_mk_median',
    'x_oct12_25_mk_median',
    'x_oct19_nov1_mk_median',
    'x_oct26_nov8_mk_median',
    'x_nov2_15_mk_median',
    'x_nov9_22_mk_median'
    };

median_names = {
    'x_oct26_nov8_mk_median',
    };

short_names = cell(size(median_names));

for i=1:length(median_names)
  median_names_i = char(median_names(i));
  short_names(i) = {median_names_i(1:end-10)};
end

for i=1:length(short_names)
  fname = ['../', bindir, '/', char(median_names(i))];
  
  % cp file to repository
  command = sprintf('cp -fp %s %s/%s', fname, repository, char(short_names(i)));
  [s,w] = my_system(command);
  

  % output slice format
  outfile = [repository, '/', char(short_names(i)), '_slice'];
  xname = char(median_names(i));
  
  cyl_ascii(outfile, xname);
  
  % output shell format
  outfile = [repository, '/', char(short_names(i)), '_shell'];
  xname = char(median_names(i));
  
  sphere_ascii(outfile, xname, r0, n_long, n_lat);
end