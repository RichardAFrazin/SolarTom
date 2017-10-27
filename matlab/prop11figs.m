%this makes the tomography figs for prop11

pth = '/Users/frazin/ML/Solar/Data/Alberto_Mar07/';

fn1 = 'x_mk4_40_60_lambda-5e-4_binfac2_2003.0805.0818';
fn2 = 'x_c2_38_90_lambda-5e-4_binfac3_2006.0609.0622';

Rmin = 1.1; Rmax = 1.85; nrad = 40; ntheta = 60; endian = 'ieee-le';
[xne1,rad1,lat1,lon1] = readtom_sph(fn1,pth,nrad,ntheta,Rmin,Rmax,endian);

Rmin = 2.3; Rmax = 6.1; nrad = 38; ntheta = 90; endian = 'ieee-le';
[xne2,rad2,lat2,lon2] = readtom_sph(fn2,pth,nrad,ntheta,Rmin,Rmax,endian);


figure(10);
subplot(2,1,1)
k=3;
im = reshape(xne1(k,:,:),60,120); 
imagesc(lon1,lat1,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['r = ',num2str(fix(100*rad1(k))/100), ' Rs    5-18 Aug 2003']);
xlabel('Carrington longitude'); ylabel('latitude');

subplot(2,1,2)
set(gca,'FontSize',18)
k = 25;
im = reshape(xne1(k,:,:),60,120); 
imagesc(lon1,lat1,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['r = ',num2str(fix(100*rad1(k))/100), ' Rs    5-18 Aug 2003']);
xlabel('Carrington longitude'); ylabel('latitude');

figure(12);
subplot(2,1,1)
k=3;
im = reshape(xne2(k,:,:),90,180); 
imagesc(lon2,lat2,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['r = ',num2str(fix(100*rad2(k))/100), ' Rs    9-22 June 2006']);
xlabel('Carrington longitude'); ylabel('latitude');

subplot(2,1,2)
k = 25;
im = reshape(xne2(k,:,:),90,180); 
imagesc(lon2,lat2,im); colorbar;
set(gca,'FontSize',18,'YTick',[-60,-30,0,30,60],'YDir','normal')
title(['r = ',num2str(fix(100*rad2(k))/100), ' Rs    9-22 June 2006']);
xlabel('Carrington longitude'); ylabel('latitude');

figpath = '/Users/frazin/docs/Tomo/Proposals/Prop11_HGI07/';
fn1 = [figpath,'fig_mk4tom.jpg'];
fn2 = [figpath,'fig_c2tom.jpg'];

str1 = ['print -f10 -djpeg30 ',fn1];
str2 = ['print -f12 -djpeg30 ',fn2];

eval(str1);
eval(str2);





