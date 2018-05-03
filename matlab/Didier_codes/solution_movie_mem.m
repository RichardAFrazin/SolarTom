function solution_movie_mem(x,rad,lon,lat);
nrad = size(rad,2);
ntheta = size(lat,2);

x=reshape(x,nrad,ntheta,2*ntheta);
figure;
colormap(jet(128));
for k = 1:nrad
  im = reshape(x(k,:,:),ntheta,2*ntheta);
  
  if (min(min(min(x))) < 0)
    imagesc(lon,lat,im); colorbar;
    title(['Ne  r(',num2str(k),') = ',num2str(rad(k))]);      
    xlabel('longitude');
    ylabel('latitude');
  else
    %imagesc(lon,lat,sqrt(im)); colorbar;
    %title(['SQRT(Ne)  r(',num2str(k),') = ',num2str(rad(k))]);
    imagesc(lon,lat,log10(im),[3 6]); colorbar;
    axis image;
    title(['LOG10(Ne)  r(',num2str(k),') = ',num2str(rad(k))]);
    xlabel('longitude');
    ylabel('latitude');
  end
  set(gca,'YDir','normal');
  pause;
end