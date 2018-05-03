function phi_movie_compare(fname1, fname2, pause_time, num_pass, logscale)
% function phi_movie_compare(fname1, fname2, pause_time, [num_pass, logscale])

if (nargin < 4)
  num_pass = 1;
  logscale  = 0;
end

[x1, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname1);
[x2, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname2);

figure(1);

set(gcf, 'DoubleBuffer', 'On');

i = 1;

if (logscale)
  x1 = log10(x1 + 1);
  x2 = log10(x2 + 1);
end

if (num_pass == 0)
  while (1)
    if (i > nphi)
      i = 1;
    end
    
    subplot(1,3,1);
    imagesc(r,z,squeeze(x1(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    drawnow;
    
    subplot(1,3,2);
    imagesc(r,z,squeeze(x2(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    
    subplot(1,3,3);
    imagesc(r,z,abs(squeeze(x1(:,i,:))' - squeeze(x2(:,i,:))'));
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
    
    subplot(1,3,1);
    imagesc(r,z,squeeze(x1(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    
    subplot(1,3,2);
    imagesc(r,z,squeeze(x2(:,i,:))');
    set(gca, 'Ydir', 'normal');
    xlabel('r');
    ylabel('z');
    title(['\phi = ', num2str(phi(i))]);
    axis('equal');
    axis('tight');
    colorbar;
    
    subplot(1,3,3);
    imagesc(r,z,abs(squeeze(x1(:,i,:))' - squeeze(x2(:,i,:))'));
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