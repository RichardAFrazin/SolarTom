function writetom_cyl(x, nrad, nphi, nz, fname)

if (prod(size(x)) ~= nrad * nphi * nz)
  error('Incorrrect dimensions');
end

fid = fopen(fname, 'w');
if (fid == -1)
  error(['Could not open file: ', fname]);
end

count = fwrite(fid, x(:), 'float32');

if (count ~= prod(size(x)))
  error(['Did not write enough info to file: ', fname]);
end

fclose(fid);