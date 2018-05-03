function i = load_fitstest(filename)
% i = function load_fitstest(filename)

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

fid = fopen(filename);

if (fid == -1)
  error(['Could not open file: ', filename]);
end

[y,r] = strtok(fliplr(filename),'x');
r = r(2:end);
y = fliplr(y);
y = sscanf(y, '%d');
x = fliplr(strtok(r,'_'));
x = sscanf(x, '%d');

[i, count] = fread(fid, x*y, 'float32');

if (count ~= x*y)
  error(['Error reading file: ', filename]);
end

fclose(fid);

i = reshape(i,y,x);