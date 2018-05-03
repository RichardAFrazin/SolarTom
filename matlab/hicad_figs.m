%this makes the tomography figs for prop11

pth = '/Users/frazin/ML/Solar/Data/Alberto_Mar07/';

%fn1 = 'x_mk4_40_60_lambda-5e-4_binfac2_2003.0805.0818';
fn1 = 'x_c2_38_90_lambda-5e-4_binfac3_2006.0609.0622';
fn2 = 'x_c2_38_90_lambda-2.5e-4_binfac3_2006.0609.0622-14files';

%Rmin = 1.1; Rmax = 1.85; nrad = 40; ntheta = 60; endian = 'ieee-le';
%[xne1,rad1,lat1,lon1] = readtom_sph(fn1,pth,nrad,ntheta,Rmin,Rmax,endian);

Rmin = 2.3; Rmax = 6.1; nrad = 38; ntheta = 90; endian = 'ieee-le';
[xne1,rad,lat,lon] = readtom_sph(fn1,pth,nrad,ntheta,Rmin,Rmax,endian);
[xne2]             = readtom_sph(fn2,pth,nrad,ntheta,Rmin,Rmax,endian);


k=3;
figure(10);
subplot(2,1,1)
im = reshape(xne1(k,:,:),90,180); 
imagesc(lon,lat,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['90 images, r = ',num2str(fix(100*rad(k))/100), ' Rs']);
xlabel('Carrington longitude'); ylabel('latitude');
v = caxis;

subplot(2,1,2)
set(gca,'FontSize',18)
im = reshape(xne2(k,:,:),90,180); 
imagesc(lon,lat,im); caxis(v); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['14 images, r = ',num2str(fix(100*rad(k))/100), ' Rs']);
xlabel('Carrington longitude'); ylabel('latitude');


figure(12);
k=33;
subplot(2,1,1)
im = reshape(xne1(k,:,:),90,180); 
imagesc(lon,lat,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['90 images, r = ',num2str(fix(100*rad(k))/100), ' Rs']);
xlabel('Carrington longitude'); ylabel('latitude');
v = caxis;

subplot(2,1,2)
set(gca,'FontSize',18)
im = reshape(xne2(k,:,:),90,180); 
imagesc(lon,lat,im); caxis(v); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['14 images, r = ',num2str(fix(100*rad(k))/100), ' Rs']);
xlabel('Carrington longitude'); ylabel('latitude');



figpath = '/Users/frazin/docs/Tomo/Papers/HiCadence/';
fnn1 = [figpath,'fig_1.jpg'];
fnn2 = [figpath,'fig_2.jpg'];

str1 = ['print -f10 -djpeg30 ',fnn1];
str2 = ['print -f12 -djpeg30 ',fnn2];
str3 = ['!cp hicad_figs.m ',figpath,'hicad_figs.m']
eval(str1);
eval(str2);
eval(str3);




