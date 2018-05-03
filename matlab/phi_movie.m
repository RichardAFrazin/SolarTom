function phi_movie(filename, pause_time, num_pass, logscale)
% function phi_movie(filename, pause_time, [num_pass, logscale])

if (nargin < 3)
  num_pass = 1;
  logscale  = 0;
end

if (nargin == 3)
  logscale = 0;
end

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename);

figure(1);

set(gcf, 'DoubleBuffer', 'On');

i = 1;

if (logscale)
  x = log10(x + 1);
end

if (num_pass == 0)
  while (1)
    if (i > nphi)
      i = 1;
    end
    
    imagesc(r,z,squeeze(x(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    drawnow;
    
    i = i + 1;
    pause(pause_time);
  end
else
  pass = 1;
  while (pass <= num_pass)
    if (i > nphi)
      i = 1;
      pass = pass + 1;
    end
    
    imagesc(r,z,squeeze(x(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    drawnow;
    
    i = i + 1;
    pause(pause_time);
  end
end