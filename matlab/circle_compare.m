function c = circle_compare(raw_list, raw_dir, xname, r0, path, ...
                            scale_median)
% function circle_compare(raw_list, raw_dir, xname, r0, [path,
%                         scale_median])

N = 256;

load('~/src/srt/eit_data/cmap.mat');
load('~/src/srt/eit_data/cmap_data.mat');

if (nargin < 6)
  scale_median = 1;
end

for i=1:length(raw_list)
  % ---------------------------------------------------------------------------
  % Compute and display image
  
  raw_fname = [raw_dir, '/', char(raw_list(i))];
  
  fid = fopen(raw_fname, 'r');
  
  if (fid == -1)
    error(['Error opening file: ', raw_fname]);
  end
  
  % header info stored as doubles at the begining of the file
  [eit_header, count] = fread(fid, 4, 'float64');
  
  if (count ~= 4)
    disp(fits_name);
    error(['Error prcoessing file header: ', raw_fname]);
  end

  center_x = eit_header(1);
  center_y = eit_header(2);
  solar_r = eit_header(3);
  exp_time = eit_header(4);
  
  [eit_image, count] = fread(fid, 'int16');

  if (count ~= 1024^2 & count ~= 512^2)
    disp(fits_name);
    error(['Error prcoessing file image: ', raw_fname]);
  end
  
  fclose(fid);
  
  im_size = sqrt(count);
  
  eit_image = reshape(eit_image, im_size, im_size);
  
  alpha = center_x - im_size/2;
  beta = center_y - im_size/2;
  
  x = linspace(-im_size/2 + alpha, im_size/2 + alpha, im_size) / solar_r;
  y = linspace(-im_size/2 + beta, im_size/2 + beta, im_size) / solar_r;

  Rmax_eit = max(max(x),max(y));

  % compute brightness scale factor
  gamma = 12.6 / exp_time * im_size / 1024;
  
  %subplot(1,2,1);
  figure(1);
  I = find(eit_image(:) <= eps);
  eit_image(I) = 10*eps;
  if (scale_median)
    imagesc(x,y,log(eit_image/median(eit_image(:))),log([7.2e2 ...
                        1.8e3]/median(eit_image(:))));
    set(gca, 'Ydir', 'normal');
  else
    imagesc(x,y,log(eit_image),log([7.2e2 1.8e3]));
    set(gca, 'Ydir', 'normal');
  end
  axis('equal');
  axis('tight');
  colormap(cmap);
  xlabel('x (Rsun)');
  ylabel('y (Rsun)');
  
  center = ([center_x, center_y] - im_size/2)/solar_r;
  NOP = 2^10;
  style = 'r-';
  hold on;
  c = circle(center, r0, NOP, style);
  set(c, 'LineWidth', 2);
  hold off;
  
  
  % ---------------------------------------------------------------------------
  % Compute (if not cached in cp_dir) and disk projection
  
  raw_date = sscanf(char(raw_list(i)), 'efz%4d%2d%2d.%2d%2d');
  
  year = raw_date(1);
  month = raw_date(2); 
  day = raw_date(3); 
  hour = raw_date(4); 
  minute = raw_date(5);
  
  command = ['~/src/srt/do_get_orbit ', ...
             num2str(year), ' ', ...
             num2str(month), ' ', ...
             num2str(day), ' ', ...
             num2str(hour), ' ', ...
             num2str(minute)];
  
  [s, stdout] = my_system(command);
    
  orbit_info = sscanf(stdout, '%g\t%g\t%g\t%g');
    
  if (length(orbit_info) ~= 4)
    error('Maformed data retuned by do_get_orbit');
  end
    
  carlong = orbit_info(1);
  sun_ob = orbit_info(2:end);
    
  if (nargin >= 5)
    [circ_proj, y, z] = circle_proj(xname, r0, N, carlong, sun_ob, path);
  else
    [circ_proj, y, z] = circle_proj(xname, r0, N, carlong, sun_ob);
  end

  delta_y = y(end) - y(end-1);
 
  extra_pix = ceil(2 * Rmax_eit / delta_y - N);
  S = N + extra_pix;
  
  circ_proj2 = zeros(S);
  n = ceil(extra_pix/2);
  
  circ_proj2(n+1:N+n,n+1:N+n) = circ_proj;
  
  p = (1:n) * delta_y + r0;
  y2 = [-fliplr(p), y, p];
  z2 = y2;

  figure(2)
  imagesc(y2, z2, sqrt(circ_proj2), log10(clim));
  set(gca, 'Ydir', 'normal');
  axis('square');
  xlabel('x (Rsun)');
  ylabel('y (Rsun)');
  colormap(cmap_data);
  
  drawnow;
end