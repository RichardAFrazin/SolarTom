function [carlong, sun_ob] = do_get_orbit(year, month, day, hour, ...
                                          minute, second)
%[carlong, sun_ob] = do_get_orbit(year, month, day, hour, minute, [second])


if (nargin == 5)
  second = 0;
end

pd = pwd;
cd('~/src/srt/');

command = ['./do_get_orbit ', ...
           num2str(year), ' ', ...
           num2str(month), ' ', ...
           num2str(day), ' ', ...
           num2str(hour), ' ', ...
           num2str(minute), ' ', ...
           num2str(second)];
    
[s, stdout] = my_system(command);

cd(pd);

orbit_info = sscanf(stdout, '%g\t%g\t%g\t%g');
    
if (length(orbit_info) ~= 4)
  error('Maformed data retuned by do_get_orbit');
end
    
carlong = orbit_info(1);
sun_ob = orbit_info(2:end);