function rad_movie(filename, pause_time, logscale)
% function rad_movie(filename, pause_time, [logscale])

if (nargin < 3)
  logscale  = 0;
end

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename);

figure(1);

set(gcf, 'DoubleBuffer', 'On');

i = 1;

if (logscale)
  x = log10(x + 1);
end

cmin = min(min(min(x)));
cmax = max(max(max(x)));
clim = [cmin cmax];

while (1)
  if (i > nrad)
    i = 1;
  end
  
  imagesc(phi,z,squeeze(x(i,:,:))',clim);
  set(gca, 'Ydir', 'normal');
  xlabel('\phi');
  ylabel('z');
  title(['r = ', num2str(r(i))]);
  axis square;
  colorbar;
  drawnow;
  
  i = i + 1;
  pause(pause_time);
end