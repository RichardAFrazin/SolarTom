function write_ascii(fname_in, fname_out)
% function write_ascii(fname_in, fname_out)

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname_in);

fid = fopen(fname_out, 'w');
if (fid == -1)
  error(['Error opening file: ', fname_out]);
end

if(0)
fprintf(fid, '# ASCII conversion of : %s\n', fname_in);
fprintf(fid, '#\n');
fprintf(fid, '# nrad = %d\tnphi = %d\tnz = %d\n', nrad, nphi, nz);
fprintf(fid, '# r_min = %f\tr_max = %f\n', r(1), r(end));
fprintf(fid, '# phi_min = %f\tphi_max = %f\n', phi(1), phi(end));
fprintf(fid, '# z_min = %f\tz_max = %f\n', z(1), z(end));
fprintf(fid, '# \n');
fprintf(fid, ['# Below are nz lines terminated by a \\n.  Each line ', ...
              ' contains nrad * nphi\n']);
fprintf(fid, ['# floating point numbers separated ', ...
              'by tabs.  Each line represents a z slice of\n']); 
fprintf(fid, ['# the ', ...
              'data cylinder.  Each line is obtained by stacking ', ...
              'a nrad x xphi matrix\n']);
fprintf(fid, ['# column by column and then ', ...
              'reading off the values in the resultant vector.  ', ...
              'The\n']);
fprintf(fid, ['# first z-slice corresponds to z = zmin.\n']);
end

for i=1:nz
  x_i = x(:,:,i);
  x_i = x_i(:);
  
  if (length(x_i) ~= nrad * nphi)
    error('Assertion failed');
  end
  
  for j=1:length(x_i)
    if (j == length(x_i))
      fprintf(fid, '%g', x_i(j));
    else
      fprintf(fid, '%g\t', x_i(j));
    end
  end
  
  fprintf(fid, '\n');
end

fclose(fid);