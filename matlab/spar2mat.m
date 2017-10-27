function S = spar2mat(fname)

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

bindir = ['../', bindir];

fid_n = fopen([bindir, 'n', fname], 'r');
n = fread(fid_n, 'int32');
fclose(fid_n);

fid_i = fopen([bindir, 'i', fname], 'r');
i = fread(fid_i, 'int32');
fclose(fid_i);

fid_v = fopen([bindir, 'v', fname], 'r');
v = fread(fid_v, 'float32');
fclose(fid_v);

i = i + 1;  % first index is now 1;

j = zeros(n(end),1);

nd = diff(n);

count = 1;

for k=1:length(nd);
  j(count:(count+nd(k))-1) = k;
  count = count + nd(k);
end

n_col = length(n) - 1;
n_row = max(max(i));

S = sparse(i,j,v,n_row,n_col);