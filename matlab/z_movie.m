function z_movie(filename, pause_time, image_size, logscale)
% function z_movie(filename, pause_time, image_size, [logscale])

if (nargin < 3)
  logscale  = 0;
end

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename);

figure(1);

set(gcf, 'DoubleBuffer', 'On');

if (logscale)
  x = log10(x + 1);
end

cmin = min(min(min(x)));
cmax = max(max(max(x)));
clim = [cmin cmax];

xcart = zeros(nz,image_size,image_size);

for i=1:nz
  i/nz
  xcart(i,:,:) = z_slice(x,i,image_size);
end

while (1)
  for i=1:nz
    xcart_i = squeeze(xcart(i,:,:));
    
    imagesc(r,r,xcart_i,clim);
    set(gca, 'Ydir', 'normal');
    title(['z = ', num2str(z(i))]);
    axis square;
    colorbar;
    drawnow;
    
    pause(pause_time);
  end
  
  for i=(nz-1):-1:2
    xcart_i = squeeze(xcart(i,:,:));
    
    imagesc(r,r,xcart_i,clim);
    set(gca, 'Ydir', 'normal');
    title(['z = ', num2str(z(i))]);
    axis square;
    colorbar;
    drawnow;
    
    pause(pause_time);
  end
end