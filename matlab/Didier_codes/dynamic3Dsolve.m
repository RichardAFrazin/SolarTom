%CVRUN 
global current_x

% path to data
path='./RFrazin_codes/bindata/';
baseproj='c2.14im';
basereg='hlaplac_nr60nt60';
objectsize=432000;

[A,y,block]  = dynamicsparsebuildproj(path,baseproj,objectsize);
[Rs,Rt]  = dynamicsparsebuildreg(path,baseproj,basereg);

% split the matrices between data set & validation set
% %choose image number to compute validation on: m
% m = 10;
% rowstart = block(m);
% rowend = block(m+1);
% Ad = [A(1:rowstart-1,:);A(rowend:size(A,1),:)];
% yd = [y(1:rowstart-1) ; y(rowend:size(A,1))];
% As = A(rowstart:rowend-1,:);
% ys = y(rowstart:rowend-1,:);
% clear A

% random pixels as validation set 
npix=size(A,1);
ls = randperm(npix);
%ls = unique(randi(npix,2*npix/10,1))-1;
ls = ls(1:npix/10);
As = A(ls,:); 
ys = y(ls,:);
Ad = A;
Ad(ls,:) = [];
yd = y;
yd(ls,:) = [];
%clear A;

% unregularised sol
[x,flag,relres,iter,resvec,relvec] = lsqr(A,y,[],200);
solution_fname='xdyn_c2.14im_noreg';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);


%lambda=1e-5; mu=1e-5;
lambda = 1e-6; mu = 1e-6;
% 1st run to find an initial solution
current_x = solve_reg_lsqr(Ad,yd,Rs,Rt,lambda,mu,50);

%optimize for the regularization parameters
%score = @(q) cvcalc(q,Ad,yd,As,ys,Rs,Rt);
score = @(q) gcvcalc(q,A,y,Rs,Rt);
regparm = [lambda,mu];
regparm_optim = fminsearch(score,regparm,...
    optimset('TolX',.3,'Display','iter',...
    'PlotFcns','optimplotx'));
lambda = regparm_optim(1);
mu = regparm_optim(2);

  
%final run with best regul params
[x,flag,relres,iter,resvec] = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,200);




%write solution
%solution_fname='xdyn_c2.14iml1e-6m1.6e-4';
%solution_fname='xdyn_c2.14iml6.5e-7m1.2e-5';
%solution_fname='xdyn_c2.14iml6.5e-7m1.2e-6';
%solution_fname='xdyn_c2.14iml6e-7m1.2e-7';
solution_fname='xdyn_c2.14iml1.9e-7m2.2e-8';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%display solution
[x,rad,lat,lon,time] = dynamic_solution_movie(solution_fname,'c2.14im',path,60,60,2.3,8);

%read projected solution
ima_fname='comp_xdyn_c2.14iml6.5e-7m1.2e-5_C2-PB-20060608_0056.Mars.dat';
fid=fopen([path,'Compare/',ima_fname],'rb');
ima1=fread(fid,[512,512],'float32');
fclose(fid);
ima_fname='orig_C2-PB-20060608_0056.Mars.dat';
fid=fopen([path,'Compare/',ima_fname],'rb');
imao=fread(fid,[512,512],'float32');
fclose(fid);
ima1(find(ima1 <0))=0;
imao(find(imao <0))=0;
cmax = max(max(max(ima1)));
subplot(1,2,1);
imagesc(realsqrt(ima1)/realsqrt(cmax)); caxis([0, 1]); set(gca,'YDir','normal');
axis square;
subplot(1,2,2);
imagesc(realsqrt(imao)/realsqrt(cmax)); caxis([0, 1]); set(gca,'YDir','normal');
axis square;

imagesc(ima1-imao); set(gca,'YDir','normal');
axis square;

%%%%%%%%%
% non neg lsqr solutions

solution_fname='xdyn_c2.14iml1e-5m1e-4';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 1e-5; mu = 1e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);


solution_fname='xdyn_c2.14iml1e-6m1.6e-4';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 1e-6; mu = 1.6e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);


solution_fname='xdyn_c2.14iml6.5e-7m1.2e-5';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 6.5e-7; mu = 1.2e-5;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

solution_fname='xdyn_c2.14iml6.5e-7m1.2e-6';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 6.5e-7; mu = 1.2e-6;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

solution_fname='xdyn_c2.14iml6e-7m1.2e-7';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 6e-7; mu = 1.2e-7;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);


solution_fname='xdyn_c2.14iml1.9e-7m2.2e-8';
[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,60,60,2.3,8.);
xini = x(:);
opts.tol=1e-4;
lambda = 1.9e-7; mu = 1.2e-8;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname=[solution_fname,'nn']
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%%%%%%%%%
%TV
clear opts;
opts.tol=1e-4;
epsilon = .02*mean(y)*sqrt(numel(y)); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
ntime=14;
xini = A'*y;
An = 1e-5; % norm of A
alpha=1; beta=1;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An);
solution_fname='xdyn_c2.14im_TVnn'
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

alpha=1; beta=10;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An);
solution_fname='xdyn_c2.14im_TVnn_1_10'
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

alpha=1; beta=100;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An);
solution_fname='xdyn_c2.14im_TVnn_1_100'
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

alpha=10; beta=1;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An);
solution_fname='xdyn_c2.14im_TVnn_10_1'
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 77 images
%%%%%%%%%%%%%%%%%%%%%%%%%%%
baseproj='c2.77im';
basereg='hlaplac_nr60nt60';
objectsize=432000;
path='./RFrazin_codes/bindata/';

[A,y,block]  = dynamicsparsebuildproj(path,baseproj,objectsize);
[Rs,Rt]  = dynamicsparsebuildreg(path,baseproj,basereg);


%normA=1e-5;
%normRs=normest(Rs,1e-3) %  ~4
%normRt=normest(Rt,1e-3) %  ~12.7
%lambda = normA/normRs; %2.5e-6
%mu = normA/normRt;     %7.86e-7
lambda=2.5e-6; mu=1e-6;
x = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,200);
solution_fname='xdyn_c2.77im_l2.5e-6_m1e-6';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname='xdyn_c2.77im_l2.5e-6_m1e-6nn';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

lambda=2.5e-6; mu=1e-5
x = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,200);
solution_fname='xdyn_c2.77im_l2.5e-6_m1e-5';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
solution_fname='xdyn_c2.77im_l2.5e-6_m1e-5nn';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

clear opts;
opts.tol=1e-4;
epsilon = .02*mean(y)*sqrt(numel(y)); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
ntime=77;
xini = A'*y;
An = 1e-5; % norm of A
alpha=1; beta=1;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An,Rt);
solution_fname='xdyn_c2.77im_TVnn';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

alpha=1; beta=10;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An,Rt);
solution_fname='xdyn_c2.77im_TVnn_1_10';
fid=fopen([path,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

[x,rad,lat,lon,time] = dynamic_solution_movie(solution_fname,baseproj,path,60,60,2.3,8);
