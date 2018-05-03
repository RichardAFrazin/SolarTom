function fix_r(fname, new_fname)
  
[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname);

x_new = zeros(size(x));

for i=1:nz
  x_new(:, :, i) = x(:, :, nz - i + 1);
end

fid = my_fopen(new_fname, 'w');
i = fwrite(fid, x_new, 'float32');
if (i ~= prod(size(x_new)))
  error('error writing to file');
end
fclose(fid);