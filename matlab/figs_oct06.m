%this makes figures from the preliminary EUV tomography results.
%  the readtom_sph routine is called with the little endian flag

%location of binary files
pth = '/Users/frazin/ML/Solar/Data/Alberto_Oct06/';
nrad = 13; ntheta = 30; Rmax = 1.26; endian = 'ieee-le';


%%Load Reconstructions

fig_offset = 0;%changes figure numbers
%[x171,rad,lat,lon] = readtom_sph('x_171_13_30_lambda-050',pth,nrad,ntheta,Rmax,endian);
%[x195] = readtom_sph('x_195_13_30_lambda-050',pth,nrad,ntheta,Rmax,endian);
%[x284] = readtom_sph('x_284_13_30_lambda-050',pth,nrad,ntheta,Rmax,endian);
%[dem] = readtom_sph('DEM_13_30_lambda-050',pth,nrad,ntheta,Rmax,endian);
%[Te] = readtom_sph('TE_13_30_lambda-030',pth,nrad,ntheta,Rmax,endian);

%fig_offset = 1000;%changes figure numbers
[x171,rad,lat,lon] = readtom_sph('x_171_13_30_lambda-600_binfac1',pth,nrad,ntheta,Rmax,endian);
[x195] = readtom_sph('x_195_13_30_8hs-bins_lambda-600_binfac1',pth,nrad,ntheta,Rmax,endian);
[x284] = readtom_sph('x_284_13_30_lambda-600_binfac1',pth,nrad,ntheta,Rmax,endian);
[dem] = readtom_sph('DEM_13_30_lambda-600_binfac1-nonzero',pth,nrad,ntheta,Rmax,endian);
[Te] = readtom_sph('TE_13_30_lambda-600_binfac1-nonzero',pth,nrad,ntheta,Rmax,endian);



ll = find(x284 < 0); x284(ll) = 0;

figure(601 + fig_offset); 
ll = reshape(Te(1,:,:),30,60); 
imagesc(lon,lat,ll); colorbar;
set(gca,'YDir','normal','Fontsize',18,'YTick',[-80,-40,0,40,80]);
xlabel('longitude (deg)');
ylabel('latitude (deg)');
title('T (10^6 K), 1.0 < r < 1.02 Rs')

figure(602 + fig_offset);
ll = reshape(Te(5,:,:),30,60);
imagesc(lon,lat,ll); colorbar;
set(gca,'YDir','normal','Fontsize',18,'YTick',[-80,-40,0,40,80]);
xlabel('longitude (deg)');
ylabel('latitude (deg)');
title('T (10^6 K), 1.1 < r < 1.12 Rs');

figure(603 + fig_offset);
subplot(1,3,1)
k = 8;
ll = reshape(Te(:,:,k),13,30);
imagesc(rad,lat,ll'); title(['\phi = ',num2str(6*(k-.5)),' deg']);
set(gca,'YDir','normal','Fontsize',18,'XTick',[1.1,1.2]);
ylabel('Latitude (deg)');  xlabel('Radius (Rs)');


subplot(1,3,2)
k = 25;
ll = reshape(Te(:,:,k),13,30);
imagesc(rad,lat,ll'); title(['\phi = ',num2str(6*(k-.5)),' deg']);
set(gca,'YDir','normal','Fontsize',18,'XTick',[1.1,1.2]);


subplot(1,3,3)
k = 53;
ll = reshape(Te(:,:,k),13,30); 
imagesc(rad,lat,ll'); title(['\phi = ',num2str(6*(k-.5)),' deg']);
set(gca,'YDir','normal','Fontsize',18,'XTick',[1.1,1.2]);
colorbar;

if 1
  print -f601 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/temp0.eps 
  print -f602 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/temp1.eps 
  print -f603 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/temp3.eps 
end


figure(201 + fig_offset); subplot(2,1,1)
ll = reshape(x171(1,:,:),30,60); 
maxll = max(max(ll));  ll = ll/maxll;
imagesc(lon,lat,ll); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('17.1 nm band emissivity, 1.0 < r < 1.02 Rs')

subplot(2,1,2);
kk = reshape(x171(5,:,:),30,60); kk = kk/maxll;
imagesc(lon,lat,kk); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('17.1 nm band emissivity, 1.1 < r < 1.12 Rs')


figure(211 + fig_offset); subplot(2,1,1)
ll = reshape(x195(1,:,:),30,60); 
maxll = max(max(ll));  ll = ll/maxll;
imagesc(lon,lat,ll); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('19.5 nm band emissivity, 1.0 < r < 1.02 Rs')

subplot(2,1,2);
kk = reshape(x195(5,:,:),30,60); kk = kk/maxll;
imagesc(lon,lat,kk); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('19.5 nm band emissivity, 1.1 < r < 1.12 Rs')


figure(221 + fig_offset); subplot(2,1,1)
ll = reshape(x284(1,:,:),30,60); 
maxll = max(max(ll));  ll = ll/maxll;
imagesc(lon,lat,ll); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('28.4 nm band emissivity, 1.0 < r < 1.02 Rs')

subplot(2,1,2);
kk = reshape(x284(5,:,:),30,60); kk = kk/maxll;
imagesc(lon,lat,kk); colorbar;
set(gca,'YDir','normal','Fontsize',18);
xlabel('Carrington longitude (deg)');
ylabel('latitude (deg)');
title('28.4 nm band emissivity, 1.1 < r < 1.12 Rs')

if 1
  print -f201 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/latlong171.eps
  print -f211 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/latlong195.eps 
  print -f221 -depsc2 /Users/frazin/docs/Tomo/Posters/AGUf06/latlong284.eps 
end

return;


%%% Filter Stuff

%read file with the EIT bandpasses
fid = fopen([pth,'eit_bandpasses_nw186.binary'],'rb');
crap = fread(fid,4*186,'float32',endian);  fclose(fid);
lambda = crap(1:186);
b171 = crap(187:372); b195 = crap(373:558); b284 = crap(559:744);
b171 = b171/max(b171); b195 = b195/max(b195); b284 = b284/max(b284);


%read file with temperature sensitivities
fid = fopen([pth,'eit_dem_kernels.binary'],'rb');
crap = fread(fid,800,'float32',endian); fclose(fid);
logT = crap(1:200);
q171 = crap(201:400); q195 = crap(401:600); q284 = crap(601:800);
q171 = q171/max(q171); q195 = q195/max(q195); q284 = q284/max(q284);

TT = 10.*ones(size(logT));
TT = TT.^logT;


figure(101); %subplot(2,1,1);
plot(lambda,b171,'b-','LineWidth',6)
axis([165 185 0 1.02]);
set(gca,'FontSize',24,'XTick',[165,175,185],'YTick',[0,.5,1])
xlabel('\lambda (nm)'); ylabel('Transmission'); 
title('Filter Transmission Function')

figure(102); %subplot(2,1,2);
plot(logT,q171,'b-','LineWidth',6);
axis([5.4 6.6 0 1.02]); 
set(gca,'FontSize',24,'Xtick',[5.5,6,6.5],'YTick',[0,.5,1]);
xlabel('log_1_0(T) (deg K)'); ylabel('Response');
title('Temperature Response of Filter');


figure(111); %subplot(2,1,1);
plot(lambda,b195,'g-','LineWidth',6)
axis([175 210 0 1.02]);
set(gca,'FontSize',24,'XTick',[180,190,200,210],'YTick',[0,.5,1])
xlabel('\lambda (nm)'); ylabel('Transmission'); 
title('Filter Transmission Function')

figure(112); %subplot(2,1,2);
plot(logT,q195,'g-','LineWidth',6);
axis([5.4 6.6 0 1.02]);
set(gca,'FontSize',24,'Xtick',[5.5,6,6.5],'YTick',[0,.5,1]);
xlabel('log_1_0(T) (deg K)'); ylabel('Response');
title('Temperature Response of Filter');


figure(121); %subplot(2,1,1);
plot(lambda,b284,'r-','LineWidth',6)
axis([260 290 0 1.02]); 
set(gca,'FontSize',24,'XTick',[260,270,280,290],'YTick',[0,.5,1])
xlabel('\lambda (nm)'); ylabel('Transmission'); 
title('Filter Transmission Function')

figure(122); %subplot(2,1,2);
plot(logT,q284,'r-','LineWidth',6);
axis([5.4 6.6 0 1.02]);
set(gca,'FontSize',24,'Xtick',[5.5,6,6.5],'YTick',[0,.5,1]);
xlabel('log_1_0(T) (deg K)'); ylabel('Response');
title('Temperature Response of Filter');


if 0
  print -f101 -depsc2 /Users/frazin/docs/Job_Applications/bandfig171.eps 
  print -f102 -depsc2 /Users/frazin/docs/Job_Applications/tempfig171.eps 
  print -f111 -depsc2 /Users/frazin/docs/Job_Applications/bandfig195.eps 
  print -f112 -depsc2 /Users/frazin/docs/Job_Applications/tempfig195.eps 
  print -f121 -depsc2 /Users/frazin/docs/Job_Applications/bandfig284.eps 
  print -f122 -depsc2 /Users/frazin/docs/Job_Applications/tempfig284.eps 
end


crap = 0; clear crap;
