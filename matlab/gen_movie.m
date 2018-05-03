function gen_movie(raw_dir, frame_dir)
% function gen_movie(raw_dir, frame_dir)
%
% raw_dir =   '/home/butala/src/srt/eit_data/fits_raw';
% frame_dir = '/home/butala/src/srt/eit_data/eit_frames';


combine_images = 1;
add_date = 1;
generate_movie = 1;


file_list = dir([frame_dir, '/', 'eit_frame_*.eps']);
eit_list = {file_list.name};

file_list = dir([frame_dir, '/', 'data_frame_*.eps']);
data_list = {file_list.name};

file_list = dir([raw_dir, '/*.raw']);
raw_list = {file_list.name};


if (length(eit_list) ~= length(data_list))
  error('# eit frames must = # data frames');
end

N = length(eit_list);



% combine images
if (combine_images)
  disp('Combining images');
  
  horz_border = 'border_32x395.png';
  vert_border = 'border_856x53.png';
  
  for i=1:N
    eit_i = char(eit_list(i));
    data_i = char(data_list(i));
  
    command = sprintf(['convert %s/%s %s/%s %s/%s +append ' ...
                       '%s/combined_%.3d.png'], ...
                      frame_dir, eit_i, ...
                      [frame_dir, '/border'], horz_border, ...
                      frame_dir, data_i, ...
                      frame_dir, i);
                      
    s = my_system(command);
    
    command = sprintf(['convert %s/%s %s/combined_%.3d.png '...
                       '-append %s/frame_%.3d.png'], ...
                      [frame_dir, '/border'], vert_border, ...
                      frame_dir, i, ...
                      frame_dir, i);
    
    s = my_system(command);
    
    disp([num2str(round(i / N * 1e3)/10), '%']);
  end
else
  disp('NOT combining images');
end


% add date
if (add_date)
  disp('Adding date');
  
  for i=1:N
    raw_date = sscanf(char(raw_list(i)), 'efz%4d%2d%2d.%2d%2d%2d');
    
    year = raw_date(1);
    month = raw_date(2); 
    day = raw_date(3); 
    hour = raw_date(4); 
    minute = raw_date(5);
    second = raw_date(6);
    
    raw_date_str = datestr(datenum(year, month, day, hour, minute, second), 31);
    
    command = sprintf(['convert %s/frame_%.3d.png ' ...
                       '-gravity North -pointsize 24 ', ...
                       '-draw "gravity north text 0,8 ''%s UTC''" ', ...
                       '-type TrueColor -depth 8 %s/dframe_%.3d.png'], ...
                      frame_dir, i, raw_date_str, frame_dir, i);
    
    s = my_system(command);
    
    disp([num2str(round(i / N * 1e3)/10), '%']);
  end
else
  disp('NOT adding date');
end


% Create movie
if (generate_movie)
  disp('Generating image_list.txt');
  
  fid = my_fopen(sprintf('%s/%s', frame_dir, 'image_list.txt'), 'w');
  
  for i=1:N
    dframe_i = sprintf('dframe_%.3d.png', i);
    
    fprintf(fid, '%s\n', dframe_i);
  end
  
  fclose(fid);
  
  
  disp('Generating movie');
  
  % create movie
  
  % Currently outpuing a motion JPEG
  
  % -z flipud
  % -g input frame size
  % -k flip r and b channels
  % -f frames per second
  
  %codec = 'mjpeg';
  %codec = 'mpeg1';
  
  %command = sprintf(['transcode --use_rgb ', ...
  %                   '-i image_list.txt -x imlist=rgb,null ', ...
  %                   '-g 944x496 -z -k ', ...
  %                   '-y ffmpeg,null -F %s -o eit_cp -f 4'], codec);
  
  command = sprintf(['transcode --use_rgb ', ...
                     '-i image_list.txt -x imlist=rgb,null ', ...
                     '-g 856x448 ', ...
                     '-y divx5,null -f 4 -o eit_cp.avi']);
  disp(command)
  
  [s, w] = my_system(command, frame_dir);
else
  disp('NOT generating movie');
end