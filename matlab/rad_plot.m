function p = rad_plot(filename, r_index, logscale)
% function rad_plot(filename, r_index, logscale)

if (nargin < 2)
  logscale  = 0;
end

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename);

figure(1);

i = r_index;

if (logscale)
  x = log10(x + 1);
end

x_i = squeeze(x(i,:,:))';
x_i = sqrt(x_i);

p = imagesc(phi,z,x_i);
set(gca, 'Ydir', 'normal');
xlabel('\phi');
ylabel('z');
title(['r = ', num2str(r(i))]);