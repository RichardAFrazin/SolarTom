function cyl_ascii(outfile, xname)
  
[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(xname);

fid = my_fopen(outfile, 'w');

for i=1:nz
  slice = x(:, :, i);
  slice = slice(:);
  
  for j=1:(nrad*nphi)
    fprintf(fid, '%g', slice(j));
    
    if (j~=nrad*nphi)
      fprintf(fid, '\t');
    end
  end
  
  fprintf(fid, '\n');
end

fclose(fid);