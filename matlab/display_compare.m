function display_compare()

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

l_in = dir(['../',bindir,'t_in_*']);
l_out = dir(['../',bindir,'t_out_*']);

if(size(l_in) ~= size(l_out))
  error(['Number of t_in files should equal the number of t_out ' ...
         'files']);
end

figure(1);
set(gcf, 'DoubleBuffer', 'On');

for i=1:length(l_in)
  fname_in = ['../', bindir, l_in(i).name ];
  fid_in = fopen(fname_in,'r');
  if(fid_in == -1)
    error(['Could not open file: ', fname_in]);
  end
  i_in = fread(fid_in,'float32');
  i_in = reshape(i_in, imsize, imsize);
  fclose(fid_in);
  
  fname_out = ['../', bindir, l_out(i).name ];
  fid_out = fopen(fname_out,'r');
  if(fid_out == -1)
    error(['Could not open file: ', fname_out]);
  end
  i_out = fread(fid_out,'float32');
  i_out = reshape(i_out, imsize,imsize);
  fclose(fid_out);
  
  cmin = min(min([i_in i_out]));
  cmax = max(max([i_in i_out]));
  clim = [cmin cmax];
  
  subplot(1,2,1);
  imagesc(i_in,clim);
  %colorbar;
  title(['Data: ', fliplr(strtok(fliplr(l_in(i).name),'_'))]);
  axis square;
  drawnow;
  
  subplot(1,2,2);
  imagesc(i_out,clim);
  %colorbar;
  title('Reconstruction');
  axis square;
  drawnow;

  pause;
end