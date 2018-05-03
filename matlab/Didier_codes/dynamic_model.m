directory='./RFrazin_codes/bindata/'

model_fname='xmodel_c2'
[x0,rad,lat,lon] = readtom_sph(model_fname,directory,60,60,2.3,8.);

% build rotating model
%%%%%%%%%%%%%%%%%%%%%%%
nimages=14;
x0d = repmat(x0,[1 1 1 nimages]);
for t=2:nimages
     x0d(:,:,:,t) = circshift(x0,[0 0 t-1]);
end
rotmodel_fname='xrotmodel_c2';
fid=fopen([directory,rotmodel_fname],'wb');
fwrite(fid,x0d(:),'float32');
fclose(fid);
base='c2.14im';
x = dynamic_solution_movie(rotmodel_fname,base,directory,60,60,2.3,8);

% same for 77 images (same rotation speed)
[r3,t3,p3]=meshgrid(rad,lat,[lon-360,lon]);
nimages=77;
x0d = repmat(x0,[1 1 1 nimages]);
step = 3*13/76;
for t=2:nimages
     loni = lon - (t-1)*step;
     [r3i,t3i,p3i]=meshgrid(rad,lat,loni);
     x0d(:,:,:,t) = interp3(r3,t3,p3,cat(3,x0,x0),r3i,t3i,p3i,'*linear');
end
rotmodel_fname='xrotmodel_c2_77';
fid=fopen([directory,rotmodel_fname],'wb');
fwrite(fid,x0d(:),'float32');
fclose(fid);
base='c2.77im';
x0d = dynamic_solution_movie(rotmodel_fname,base,directory,60,60,2.3,8);
clear p3 t3 r3 r3i t3i p3i loni



% add sharp variation model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tvar = nimages/2;
x0ds=x0d;
x0ds(:,:,:,tvar:nimages) = x0d(:,:,:,tvar:nimages)*.5 ; % remove 50%
sharpmodel_fname='sharpmodel_c2';
fid=fopen([directory,sharpmodel_fname],'wb');
fwrite(fid,x0ds(:),'float32');
fclose(fid);

x = dynamic_solution_movie(sharpmodel_fname,'c2.14im',directory,60,60,2.3,8);

%%% build images from dynamic model 
%%% then test static rec
%base='c2.14im';
base='c2.77im';
objectsize=432000;
[A,y,block]  = dynamicsparsebuildproj(directory,base,objectsize);
clear y;

x0d = dynamic_readtom_sph(rotmodel_fname,base,directory,60,60,2.3,8.);
x0d=x0d(:);
y = A*x0d(:);
% add noise
%yb = y +  .05*mean(y)*randn(size(y));
yb = y +  .02*mean(y)*randn(size(y));

% x0ds = dynamic_readtom_sph(sharpmodel_fname,base,directory,60,60,2.3,8.);
% x0ds=x0ds(:);
% y = A*x0ds(:);
% % add noise
% %yb = y +  .05*mean(y)*randn(size(y));
% yb = y +  .02*mean(y)*randn(size(y));
% 

% load static A
clear A
[A,y]=staticsparseread(directory,base,432000);
clear y;

Rt0 = zeros(0,size(A,2)); % empty time regul matrix
Rs0 = zeros(0,size(A,2)); % empty spatial regul matrix
[d2r,d2theta,d2phi,Lh,Lr,L] = spatialsparsebuildreg(60,60,120,2.3,8); 
clear d2r d2theta d2phi 
clear Lh Lr 

% lsqr solution
lambda =1e-5;
[x,flag,~,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
%solution_fname='xrotmodelrec_c2_l1e-5';
solution_fname='xsharpmodelrec_c2_l1e-5';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,L,Rt0,lambda,0,xini,opts,1e-5);
%solution_fname='xrotmodelrec_c2_l1e-5_nn';
solution_fname='xsharpmodelrec_c2_l1e-5_nn';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

% gcv to compute solution%     
[m,n]=size(A);
u = randi(2,m,1);
u = 2*u-3; % -1,+1 random vector
score = @(q) gcvcalc([1e-5*q,0],A,yb,L,Rt0,u);
lambdamin=1e-5 ; lambdamax=1e2;
regparm_optim = fminbnd(score,lambdamin,lambdamax,...
    optimset('Display','iter'));
lambda = regparm_optim*1e-5;

[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
%solution_fname='xrotmodelrec_c2_lgcv'; % l=6.533e-7   
solution_fname='xsharpmodelrec_c2_lgcv';% l=6.4758e-7
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,L,Rt0,lambda,0,xini,opts,1e-5);
%solution_fname='xrotmodelrec_c2_lgcv_nn';
solution_fname='xsharpmodelrec_c2_lgcv_nn';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

% TV(theta,phi) regularized solution using TFOCS
clear opts;
opts.tol=1e-4;
epsilon = .02*mean(yb)*sqrt(numel(y)); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
xini = A'*yb;
An = 1e-5; % norm of A
[x,out,optsout] = solve_TVthetaphi(A,yb,epsilon,nrad,ntheta,nphi,xini,opts,An);
%solution_fname='xrotmodelrec_c2_TV';
solution_fname='xsharpmodelrec_c2_TV';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%%% VISU
x=solution_movie(solution_fname,directory,60,60,2.3,8.);
x = readtom_sph(solution_fname,directory,60,60,2.3,8.);
solution_movie_mem(x,rad,lon,lat);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dynamic reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,y,block]  = dynamicsparsebuildproj(directory,base,objectsize);
basereg='hlaplac_nr60nt60';
[Rs,Rt]  = dynamicsparsebuildreg(directory,base,basereg);

norm_x0d = norm(x0d,'fro');
err = @(x) norm(vec(x)-vec(x0d))/norm_x0d;

irad6 = 39;
irad3 = 8;

normA=normest(A,1e-3) % ~1e-5
normRs=normest(Rs,1e-3) %  ~4
normRt=normest(Rt,1e-3) %  ~2
lambda = normA/normRs;
mu = normA/normRt;
% 1st run to find an initial solution
x = solve_reg_lsqr(A,yb,Rs,Rt,lambda,mu,50);
%solution_fname='xrotmodelrecdyn_c2_lnorm';
solution_fname='xrotmodelrecdyn_c2_77_lnorm';
%solution_fname='xsharpmodelrecdyn_c2_lnorm';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,Rs,Rt,lambda,mu,xini,opts,1e-5);
snr4D(x,x0d,irad3,irad6,nrad,ntheta,nphi,ntime) % snr=2.02 (rot)
%solution_fname='xrotmodelrecdyn_c2_lnorm_nn';
solution_fname='xrotmodelrecdyn_c2_77_lnorm_nn';
%solution_fname='xsharpmodelrecdyn_c2_lnorm_nn';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

% gcv
[m,n]=size(A);
u = randi(2,m,1);
u = 2*u-3; % -1,+1 random vector
score = @(q) gcvcalc(1e-5*q,A,yb,Rs,Rt,u);
regparm = [.25,.5];
regparm_optim = fminsearch(score,regparm,...
    optimset('TolX',.3,'Display','iter',...
    'PlotFcns','optimplotx'));
lambda = regparm_optim(1)*1e-5;
mu = regparm_optim(2)*1e-5;
[x,flag,res,iter] = solve_reg_lsqr(A,yb,Rs,Rt,lambda,mu,200);
% lambda=2e-7, mu=5.76e-6
%solution_fname='xrotmodelrecdyn_c2_gcv';
solution_fname='xrotmodelrecdyn_c2_77_gcv';
%solution_fname='xsharpmodelrecdyn_c2_gcv';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,Rs,Rt,lambda,mu,xini,opts,1e-5);
snr4D(x,x0d,irad3,irad6,nrad,ntheta,nphi,ntime) % snr=1.7 (rot)
%solution_fname='xrotmodelrecdyn_c2_gcv_nn';
solution_fname='xrotmodelrecdyn_c2_77_gcv_nn';
%solution_fname='xsharpmodelrecdyn_c2_gcv_nn';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

% TV solution
% 
clear opts;
opts.tol=1e-4;
opts.errFcn     = @(f,dual,primal) err(primal);
epsilon = .02*mean(yb)*sqrt(numel(yb)); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
%ntime=14;
ntime=77;
xini = A'*yb;
An = 1e-5; % norm of A
%alpha=1; beta=1;
%alpha=1; beta=10;
%alpha=1; beta=100;
alpha=1; beta=1e5;
[x,out,optsOut] = solve_TVthetaphi_TVtime(A,yb,epsilon,alpha,beta,nrad,ntheta,nphi,ntime,xini,opts,An);
%solution_fname='xrotmodelrecdyn_c2_TV';
%solution_fname='xrotmodelrecdyn_c2_TV_1_10';
%solution_fname='xrotmodelrecdyn_c2_TV_1_100';
%solution_fname='xrotmodelrecdyn_c2_77_TV_1_10';
%solution_fname='xrotmodelrecdyn_c2_77_TV_1_10_testing';
%solution_fname='xrotmodelrecdyn_c2_77_TV_1_100';
solution_fname='xrotmodelrecdyn_c2_77_TV_1_1e5';
snr4D(x,x0d,irad3,irad6,nrad,ntheta,nphi,ntime) 
% 1_1  -> snr=1.93
% 1_10 -> snr=1.93  | 77im snr=1.54
% 1_100 -> snr=1.93 | 77im snr=1.4
%solution_fname='xsharpmodelrec_c2_TV';
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);

[x,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,base,directory,60,60,2.3,8.);
x = dynamic_solution_movie(solution_fname,base,directory,60,60,2.3,8);

