function map_jspace(cpu_i, M)

warning off;
addpath('~/src/srt/matlab');
warning on;

N = 13;
jspace_file = ['jspace_', num2str(cpu_i), '.mat'];

if (mod(N,2) ~= 1)
  error('N must be odd');
end

x_jspace = ['x_jspace_', num2str(cpu_i)];
sys_matrix_name = 'oct5_18_mk_block_c';
norm_matrix_name = 'oct5_18_mk_block';

lambda0 = [1.13393e-06,     7.92391e-06,     2.55166e-06];

x_infile = 'x_mk0';

jspace = zeros(N,N);

% Compute initial point
if (~exist(['./bindata/', x_jspace], 'file'))
  datestr(now)
  disp('Computing initial point');
  cd('~/src/srt');
  command = sprintf('./auto_x0 -v %s %s %s %g %g %g', ...
                    x_infile, x_jspace, sys_matrix_name, ...
                    lambda0(1), lambda0(2), lambda0(3));
  
  disp(command);
  
  [s,w] = my_system(command);
else
  disp(['Initial point: ', x_jspace]);
end

% Compute costs
l1 = logspace(log10(lambda0(1)/10000), log10(lambda0(1)*1), N);
%l2 = logspace(log10(lambda0(2)/10), log10(lambda0(2)*10), N);
l2 = lambda0(2);
l3 = logspace(log10(lambda0(3)/1000), log10(lambda0(3)*10), N);

count = 1;

for i=1:length(l1)
  %    for j=1:length(l2)
  for k=1:length(l3)
    if (mod(i * length(l3) + k, M) + 1 == cpu_i)
      disp(['i = ', num2str(i), '   k = ', num2str(k)]);
      
      datestr(now)
      
      x_jspace_ijk = sprintf('%s_%d%d', x_jspace, i, k);
      
      cd('~/src/srt');
      command = sprintf('./do_cvcalc %g %g %g %s %s %s %s', ...
                        l1(i), l2, l3(k), ...
                        norm_matrix_name, sys_matrix_name, ...
                        x_jspace, x_jspace_ijk);
      disp(command);
      
      [s,w] = my_system(command);
      
      I = findstr('normval: ', w);
      
      if (length(I) ~= 1)
        error('Unexpected output returned from callsolve_fess');
      end
      
      normval = sscanf(w(I(1):end), 'normval: %g');
      
      jspace(i,k) = normval;
    end
  end
  
  disp([num2str(round(count/N^2 * 1e3)/10), '%']);
  count = count + 1;
  %   end
end

cd('~/src/srt/matlab');
save(jspace_file, 'jspace', 'l1', 'l2', 'l3');
