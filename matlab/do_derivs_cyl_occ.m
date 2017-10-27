% do_derivs_cyl_occ.m
%
% Automate the process of generating the regularization matricies.
%
% Mark Butala
% Jul. 2004

DELTA_VAL = 200;

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

bindir = ['../', bindir];

fname = [num2str(nrad), num2str(nphi), num2str(nz), '.', ...
         num2str(rmin*100), '.', num2str(rmax*100)];



% =============================================================================
% Handle regularization matricies.

if (exist([bindir,'n','d2r',fname],'file') == 2 & ...
    exist([bindir,'i','d2r',fname],'file') == 2 & ...
    exist([bindir,'v','d2r',fname],'file') == 2 & ...
    exist([bindir,'y','d2r',fname],'file') == 2 & ...
    exist([bindir,'delta_','d2r',fname],'file') == 2 & ...
    exist([bindir,'n','dphi',fname],'file') == 2 & ...
    exist([bindir,'i','dphi',fname],'file') == 2 & ...
    exist([bindir,'v','dphi',fname],'file') == 2 & ...
    exist([bindir,'y','dphi',fname],'file') == 2 & ...
    exist([bindir,'delta_','dphi',fname],'file') == 2 & ...
    exist([bindir,'n','d2z',fname],'file') == 2 & ...
    exist([bindir,'i','d2z',fname],'file') == 2 & ...
    exist([bindir,'v','d2z',fname],'file') == 2 & ...
    exist([bindir,'y','d2z',fname],'file') == 2 & ...
    exist([bindir,'delta_','d2z',fname],'file') == 2)
  disp('do_derivs_cyl_occ: files already exist...exiting without generating data');
else
  derivs_cyl_occ(nrad, nphi, nz, rmin, grid_rmax, bindir, fname);
  
  disp('no worries, automatically running row_to_col');

  mat1 = {'d2r', 'd2z', 'dphi'};
  
  rows = nrad * nphi * nz;
  cols = rows;
  
  for i=1:length(mat1)
    fname_ij = [char(mat1(i)), fname];
    
    command_string = ['../row_to_col ', ...
                      bindir, ' ', ...
                      fname_ij, ' ', ...
                      num2str(rows), ' ', ...
                      num2str(cols)];
    
    [s,w] = my_system(command_string);
  end
  
  disp(['Writing DELTA_VAL = ', num2str(DELTA_VAL), ...
        ' to all derivative matrix delta vectors']);

  delta = DELTA_VAL * ones(nrad*nphi*nz, 1);
  
  for i=1:length(mat1)
    fname_ij = [bindir, '/', 'delta_', char(mat1(i)), fname];
    
    fid = my_fopen(fname_ij, 'w');
    count = fwrite(fid, delta, 'float32');
    if (count ~= nrad*nphi*nz)
      error(['Error writing to file: ', fname_ij]);
    end
  end
end


% =============================================================================
% Handle initial guess vector.

full_x_infile = [bindir,x_infile];

if (exist(full_x_infile,'file') == 2)
  disp(['do_derivs_cyl_occ: ', full_x_infile, [' exists, not ' ...
                      'generating '] 'new initial x']);
else
  x = zeros(1,nrad * nphi * nz);
  
  fid = fopen(full_x_infile,'wb');
  
  if(fid == -1)
    error(['unable to open file for writing: ', full_x_infile]);
  end
  
  fwrite(fid, x, 'float32');
  fclose(fid);
  
  clear x;
end

full_cv_infile = [bindir,cv_infile];

if (exist(full_cv_infile,'file') == 2)
  disp(['do_derivs_cyl_occ: ', full_cv_infile, ' exists, not generating new initial cv x']);
else
  command = ['cp ', full_x_infile, ' ', full_cv_infile];
  
  my_system(command);
end

