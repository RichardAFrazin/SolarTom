function gen_frames(raw_dir, frame_dir, cp_dir)
% function gen_frames(raw_dir, frame_dir, cp_dir)
%
% raw_dir =   '/home/butala/src/srt/eit_data/fits_raw';
% frame_dir = '/home/butala/src/srt/eit_data/eit_frames';
% cp_dir =    '/home/butala/src/srt/eit_data/fits_cp';


file_list = dir([raw_dir, '/*.raw']);
raw_list = {file_list.name};

load('~/src/srt/eit_data/cmap.mat');
load('~/src/srt/eit_data/cmap_data.mat');

cmap_data = [[204 204 204]/255; cmap_data];


N = 256;
r0 = 1.2;

interval_names = {'x_oct5_18_mk_median';
                  'x_oct12_25_mk_median';
                  'x_oct19_nov1_mk_median';
                  'x_oct26_nov8_mk_median';
                  'x_nov2_15_mk_median';
                  'x_nov9_22_mk_median'};

intervals = [2003 10 5;
             2003 10 18;
             2003 10 12;
             2003 10 25;
             2003 10 19;
             2003 11 1;
             2003 10 26;
             2003 11 8;
             2003 11 2;
             2003 11 15;
             2003 11 9;
             2003 11 22];

if (mod(length(intervals),2) ~= 0)
  error('Interval list length must be even');
end

if (length(intervals) / 2 ~= length(interval_names))
  error('Size mismatch between intervals and intercal_names');
end

int_nums = datenum(intervals(:,1), intervals(:,2), intervals(:,3));

center_time = zeros(length(intervals)/2,1);

j = 1;
for i=1:2:length(intervals)
  center_time(j) = (int_nums(i+1) - int_nums(i))/2 + int_nums(i);
  
  j = j + 1;
end

figure(1);
set(gcf, 'DoubleBuffer', 'On');
set(gcf, 'InvertHardcopy', 'Off');
set(gca, 'FontSize', 16);

figure(2);
set(gcf, 'DoubleBuffer', 'On');
set(gcf, 'InvertHardcopy', 'Off');
set(gca, 'FontSize', 16);

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
  
  size = sqrt(count);
  
  eit_image = reshape(eit_image, size, size);
  
  alpha = center_x - size/2;
  beta = center_y - size/2;
  
  x = linspace(-size/2 + alpha, size/2 + alpha, size) / solar_r;
  y = linspace(-size/2 + beta, size/2 + beta, size) / solar_r;

  Rmax_eit = max(max(x),max(y));
  
  % compute brightness scale factor
  gamma = 12.6 / exp_time * size / 1024;
  
  % the above didn't work too well
  gamma = median(eit_image(:));
  
  %subplot(1,2,1);
  figure(1);
  I = find(eit_image(:) <= eps);
  eit_image(I) = 10*eps;
  imagesc(x,y,flipud(log(eit_image/gamma)),log([7.2e2 1.8e3]/gamma));
  set(gca,'YDir','normal');
  axis('equal');
  axis('tight');
  colormap(cmap);
  h = xlabel('solar $x$ (R$_\odot$)');
  set(h, 'interpreter', 'latex');
  h = ylabel('solar $y$ (R$_\odot$)');
  set(h, 'interpreter', 'latex');
  set(gca, 'FontSize', 16);
  
  center = ([center_x, center_y] - size/2)/solar_r;
  NOP = 2^10;
  style = 'r-';
  hold on;
  c = circle(center, r0, NOP, style);
  set(c, 'LineWidth', 2);
  hold off;
  
  
  % ---------------------------------------------------------------------------
  % Compute (if not cached in cp_dir) and disk projection
  
  
  
  cp_fname = [cp_dir, '/', char(raw_list(i)), ...
              '_', num2str(r0), '_cp.mat'];
  
  if (exist(cp_fname, 'file'))
    % cp cached - load it
    load(cp_fname);
  else
    raw_date = sscanf(char(raw_list(i)), 'efz%4d%2d%2d.%2d%2d');
  
    year = raw_date(1);
    month = raw_date(2); 
    day = raw_date(3); 
    hour = raw_date(4); 
    minute = raw_date(5);
    
    command = ['../do_get_orbit ', ...
               num2str(year), ' ', ...
               num2str(month), ' ', ...
               num2str(day), ' ', ...
               num2str(hour), ' ', ...
               num2str(minute)];
    
    [s, stdout] = system(command);
    
    if (s ~= 0)
      error(['Error running command: ', command]);
    end
    
    orbit_info = sscanf(stdout, '%g\t%g\t%g\t%g');
    
    if (length(orbit_info) ~= 4)
      error('Maformed data retuned by do_get_orbit');
    end
    
    carlong = orbit_info(1);
    sun_ob = orbit_info(2:end);
    
    [min_date,I] = min(abs(datenum(year, month, day, hour, minute, 0) - ...
                           center_time));
    
    circ_proj = zeros(N);
    
    [circ_proj, y, z] = circle_proj(char(interval_names(I)), r0, N, ...
                                    carlong, sun_ob);
    
    delta_y = y(end) - y(end-1);
 
    extra_pix = ceil(2 * Rmax_eit / delta_y - N);
    S = N + extra_pix;
    
    circ_proj2 = zeros(S);
    n = ceil(extra_pix/2);
    
    circ_proj2(n+1:N+n,n+1:N+n) = circ_proj;
    circ_proj = circ_proj2;
    
    save(cp_fname, 'circ_proj', 'y', 'z');
  end
  
  delta_y = y(end) - y(end-1);
 
  extra_pix = ceil(2 * Rmax_eit / delta_y - N);
  S = N + extra_pix;
  
  n = ceil(extra_pix/2);
    
  p = (1:n) * delta_y + r0;
  y2 = [-fliplr(p), y, p];
  z2 = y2;

  max_val = 2e4;
  
  for j=1:length(y2)
    for k=1:length(z2)
      if (sqrt(y2(j)^2 + z2(k)^2) > r0)
        circ_proj(k, j) = 0;
      else
        if (circ_proj(k, j) < 1e5)
          circ_proj(k, j) = 1e5;
        end
      end
    end
  end
  
  figure(2)
  imagesc(y2,z2,sqrt(circ_proj), [0, max_val]);
  set(gca,'YDir','normal');
  axis('square');
  h = xlabel('solar $x$ (R$_\odot$)');
  set(h, 'interpreter', 'latex');
  h = ylabel('solar $y$ (R$_\odot$)');
  set(h, 'interpreter', 'latex');
  colormap(cmap_data);
  set(gca, 'FontSize', 16);
  
  drawnow;

  figure(1);
  print('-depsc2', '', sprintf('%s/eit_frame_%.3d.eps', frame_dir, i));
  
  figure(2);
  print('-depsc2', '', sprintf('%s/data_frame_%.3d.eps', frame_dir, i));

  disp([num2str(round(i/length(raw_list)*1000)/10), '% complete']);
end


figure(1);
set(gcf, 'InvertHardcopy', 'On');

figure(2);
set(gcf, 'InvertHardcopy', 'On');