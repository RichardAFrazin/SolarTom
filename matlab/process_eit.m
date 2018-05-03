function process_eit(fits_dir, raw_dir)
% function process_eit(fits_dir, raw_dir)


list = dir([fits_dir, '/efz*']);
fits_list = {list.name};

temp_file = '/tmp/process_eit.tmp';

for i=1:length(fits_list)
  fits_name = char(fits_list(i));
  
  fits_date = sscanf(fits_name, 'efz%4d%2d%2d.%2d%2d%2d');
  year = fits_date(1);
  month = fits_date(2);
  day = fits_date(3);
  hour = fits_date(4);
  minute = fits_date(5);
  second = fits_date(6);
  
  dn = datenum(year, month, day, hour, minute, second);
  
  command = sprintf('~/src/srt/eit_image %s/%s %s/%s.raw', ...
                    fits_dir, fits_name, ...
                    raw_dir, fits_name);
  
  if (dn >= datenum(2003, 11, 14, 15, 00, 10))
     command = [command, ' 512'];
  end
  
  s = system(command);
  
  if (s ~= 0)
    error(['Error executing command: ', command]);
  end
end