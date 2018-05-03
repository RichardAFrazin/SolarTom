function phi_plot(filename, phi_index, logscale)
% function phi_plot(filename, phi_index, [logscale])

if (nargin < 2)
  logscale  = 0;
end

[x, r, phi, z, nrad, nphi, nz] = readtom_cyl(filename);

%figure(2);
figure;

i = phi_index;

if (logscale)
  x = log10(x + 1);
end

imagesc(r,z,squeeze(x(:,i,:))');
set(gca, 'Ydir', 'normal');
xlabel('r');
ylabel('z');
title(['\phi = ', num2str(phi(i))]);
axis square;
colorbar;