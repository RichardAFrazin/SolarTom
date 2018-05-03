function [x,rad,lat,lon,time] = dynamic_solution_movie(filename,baseproj,path,nrad,ntheta,Rmin,Rmax,endian);
%function [x,rad,lat,lon,time] = dynamic_solution_movie(filename,path,nrad,ntheta,Rmin,Rmax,endian);
%
% This calls dynamic_readtom_sph.m and makes a movie showing the solution at 
%    constant radius.
%
%x = output array (nrad,ntheta,2*ntheta)
%r = radial grid  linspace(1+(Rmax-1)/nrad,Rmax,nrad)
%lat = latitude grid in deg (ntheta)
%lon = longitude grid in deg (2*ntheta)
%filename - dynamic solution filename
%baseproj - filename base of projection matrix
%path - include final slash
%nrad = # of radial grid points
%ntheta = # of phi grid points
%Rmax = max radius of computation grid in Rsun
%endian = an optional argument for endian 
%    or other optional format string. usually 'ieee-be' or 'ieee-le'
%
% R.Frazin & D.Vibert 12/9/2011 (from solution_movie.m of RFrazin)

if (nargin == 8)
  endianstring = endian;
else 
   endianstring = 'ieee-le';
end

[x,rad,lat,lon,time] = dynamic_readtom_sph(filename,baseproj,path,nrad,ntheta,Rmin,Rmax,endianstring);
ntime=length(time);

% if (min(x(:)) < 0)
%    disp('Warning: non-positive solution.  linear scale displayed.');
% end

figure('Position',[1,1,900,400]);
colormap(jet(128));

% loop on time
% for k = 1:ntime
%     im1 = reshape(x(1,:,:,k),ntheta,2*ntheta);
%     im10 = reshape(x(10,:,:,k),ntheta,2*ntheta);
%     im20 = reshape(x(20,:,:,k),ntheta,2*ntheta);
%     im30 = reshape(x(30,:,:,k),ntheta,2*ntheta);
%     if (k == 1 )
%         clim1  = [0.,max(max(im1))];
%         clim10 = [0.,max(max(im10))];
%         clim20 = [0.,max(max(im20))];
%         clim30 = [0.,max(max(im30))];
%     end
%     clf;
%     subplot(2,2,1);  
%     imagesc(lon,lat,im1,clim1); colorbar;
%     title(['Ne  r(',num2str(1),') = ',num2str(rad(1)),'  time = ',num2str(time(k))]);       
%     subplot(2,2,2);  
%     imagesc(lon,lat,im10,clim10); colorbar;
%     title(['Ne  r(',num2str(10),') = ',num2str(rad(10))]); 
%     subplot(2,2,3);  
%     imagesc(lon,lat,im20,clim20); colorbar;
%     title(['Ne  r(',num2str(20),') = ',num2str(rad(20))]); 
%     subplot(2,2,4);  
%     imagesc(lon,lat,im30,clim30); colorbar;
%     title(['Ne  r(',num2str(30),') = ',num2str(rad(30))]);
% 
%     set(gca,'YDir','normal');
%     pause;
% end
%close;

x=max(x,zeros(size(x)));

itime = round(((1:4)-1)*((ntime-1)/3)+1);
clim=[3 6];
for k = 1:nrad
  im1 = reshape(x(k,:,:,itime(1)),ntheta,2*ntheta);
 im10 = reshape(x(k,:,:,itime(2)),ntheta,2*ntheta);
 im20 = reshape(x(k,:,:,itime(3)),ntheta,2*ntheta);
 im30 = reshape(x(k,:,:,itime(4)),ntheta,2*ntheta);

 clf;
 subplot(2,2,1);
 imagesc(lon,lat,log10(im1),clim);colorbar;
 title(['LOG10(Ne)   time(',num2str(itime(1)),') = ',num2str(time(itime(1))),...
     '  r(',num2str(k),') = ',num2str(rad(k))]);
 set(gca,'YDir','normal');axis image;
 subplot(2,2,2);
 imagesc(lon,lat,log10(im10),clim); colorbar;
 title(['LOG10(Ne)   time(',num2str(itime(2)),') = ',num2str(time(itime(2)))]);
 set(gca,'YDir','normal');axis image;
 subplot(2,2,3);
 imagesc(lon,lat,log10(im20),clim); colorbar;
 title(['LOG10(Ne)   time(',num2str(itime(3)),') = ',num2str(time(itime(3)))]);
 set(gca,'YDir','normal');axis image;
 subplot(2,2,4);
 imagesc(lon,lat,log10(im30),clim); colorbar;
 title(['LOG10(Ne)   time(',num2str(itime(4)),') = ',num2str(time(itime(4)))]);
 set(gca,'YDir','normal');axis image;
 
 pause;
end


return;
