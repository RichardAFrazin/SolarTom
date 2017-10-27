function inspect_mk4()

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;
 
bindir = ['../', bindir];

conf_list = {
    'MK4oct03_5_18',
    'MK4oct03_12_25',
    'MK4oct03_19_1',
    'MK4oct03_26_8',
    'MK4nov03_2_15',
    'MK4nov03_9_22'
    };

fits_list = {};

for i=1:length(conf_list)
  fname = [bindir, '/', char(conf_list(i))];
  fid = my_fopen(fname, 'r');
  
  n = sscanf(fgets(fid), '%d');
  
  if (length(n) ~= 1)
    error(['Error processing: ', fname]);
  end
  
  for j=1:n
    l = fgets(fid);
    
    if (l(end) ~= 's')
      l = l(1:end-1);
    end
    
    if (l == -1)
      error(['Error processing: ', fname]);
    end
    
    fits_list(end+1) = {l};
  end
  
  fclose(fid);
end

z = (2 * grid_rmax / nz / 2) * (2 * (1:nz) - 1) - grid_rmax;

figure(1);
set(gcf, 'DoubleBuffer', 'On');

for i=1:length(fits_list)
  fits_i = char(fits_list(i));
  
  command = sprintf('../fitstest -m -s %s/%s %s/%s', ...
                    '../mk4_data', fits_i, '/tmp', 'out');
  disp(command);
  
  [s,w] = my_system(command);
  
  fid = my_fopen('/tmp/out_960x960', 'r');
  x = fread(fid, 'float32');
  fclose(fid);
  
  x = reshape(x,960,960);

  i
  
  figure(1);
  imagesc(z, z, x);
  set(gca, 'Ydir', 'normal');
  cm = colormap;
  axis('equal');
  axis('tight');
  colorbar;
  drawnow;
end
