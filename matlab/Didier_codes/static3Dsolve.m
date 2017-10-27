% 3D static reconstruction with RFrazin matrices
lambda = 1e-5; % regularization param
%lambda = 6e-7;
directory='./RFrazin_codes/bindata/'
base='c2mars.pb.2006.0608.0621.rmax8.0nr60nt60';
%base='c2.14im';

[A,y]=staticsparseread(directory,base,432000);

% unregularized solution
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[x,flag,relres,iter,resvec,relvec] = lsqr(A,y,[],200);
% relres=norm(y-A*x)/norm(y)=1.4 %
% resvec %-> ~15
% norm(x(:)) %-> ~2e+8

Rt0 = zeros(0,size(A,2)); % empty time regul matrix
Rs0 = zeros(0,size(A,2)); % empty spatial regul matrix

%build positive unregularized solution using TFOCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x0 = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs0,Rt0,0,0,x0,opts,1e-5);

% build thikonov regularized solution (d2theta,d2phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[L,y0]=staticsparseread(directory,'hlaplac_nr60nt60',432000);
[d2r,d2theta,d2phi,Lh,Lr,L] = spatialsparsebuildreg(60,60,120,2.3,8); 
%L=d2phi;
clear d2r d2theta d2phi 
clear Lh Lr 

%H=[A;lambda*L];
% y=[y;y0]
% clear L
% clear y0
% x=lsqr(H,y,[],200);


% gcv to compute solution%     
[m,n]=size(A);
u = randi(2,m,1);
u = 2*u-3; % -1,+1 random vector
score = @(q) gcvcalc([1e-5*q,0],A,y,L,Rt0,u);
lambdamin=1e-5 ; lambdamax=1e2;
regparm_optim = fminbnd(score,lambdamin,lambdamax,...
    optimset('Display','iter'));
lambda = regparm_optim*1e-5;
%nA = 1e-5;A = A/nA;lambda=lambda/nA;
[x,flag,res,iter] = solve_reg_lsqr(A,y,L,Rt0,lambda,0,200);
%x = x/nA;lambda=lambda/nA;A=A*nA;
%[x,flag,res,iter] = solve_reg_lsqr(A,y,L,Rt0,lambda,0,50);

%build positive solution using TFOCS
xini = x;
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,y,L,Rt0,lambda,0,xini,opts,1e-5);

% TV(theta,phi) regularized solution using TFOCS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear opts;
opts.tol=1e-4;
%epsilon = .05*mean(y)*sqrt(254040); % sigma*sqrt(n) ~ 30
epsilon = .02*mean(y)*sqrt(1378196); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
xini = A'*y;
An = 1e-5; % norm of A
[x,out,optsout] = solve_TVthetaphi(A,y,epsilon,nrad,ntheta,nphi,xini,opts,An);
x(x<0)=0;% dirty positivity

solution_fname='xstat_c2.14im'; % without regularization
solution_fname='xstat_c2.14im_nn'; % without regularization >0
solution_fname='xstat_c2.14iml1e-5';
solution_fname='xstat_c2.14iml1e-5nn';
%solution_fname='xstat_c2.14iml7e-7'; %l=7.33e-7 gcv solution
solution_fname='xstat_c2.14iml7e-7nn';
%solution_fname='xstat_c2.14iml1e-5-ctheta';
%solution_fname='xstat_c2.14iml1e-7-ctheta'; % l=1e-7 gcv solution
%solution_fname='xstat_c2.14iml5e-5-ctheta-only';  
%solution_fname='xstat_c2.14iml4e-6-ctheta-only'; %l=4e-6 gcv solution
%solution_fname='xstat_c2.14iml1e-5-3dlaplac';
solution_fname='xstat_c2.14iml1e-5-temp';
solution_fname='xstat_c2.14imTV-temp';

% 77 images
solution_fname='xstat_c2.77iml1e-5';
solution_fname='xstat_c2.77iml1e-5nn';
solution_fname='xstat_c2.77iml1.4e-6'; %l=1.36e-6 gcv solution
solution_fname='xstat_c2.77iml1.4e-6nn'; 
solution_fname='xstat_c2.77imTV';
solution_fname='xstat_c2.77imTVnn';

fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);
solution2ascii(directory,solution_fname,60,60,2.3,8.,6.);

[x,rad,lat,lon] = readtom_sph(solution_fname,directory,60,60,2.3,8.);
[x,rad,lat,lon]=solution_movie(solution_fname,directory,60,60,2.3,8.);

solution_movie_mem(max(x,zeros(size(x))),rad,lon,lat);
solution_movie_mem(x,rad,lon,lat);

x=reshape(x,60,60,120);

[x,rad,lat,lon]=solution_movie('x_c2.14iml1.e-5','./RFrazinMatrices/',60,60,2.3,8.);

% 3D dynamic reconstruction with RFrazin matrices
lambda = 1e-5; % spatial regularization param
mu = lambda; % temporal regularization param
[H,y]  = dynamicsparsebuildproj('./RFrazinMatrices/','c2.14im',432000);

[L,D] = dynamicsparsebuildreg('./RFrazinMatrices/','c2.14im','hlaplac_nr60nt60');
H=[H;lambda*L;mu*D];
y=[y;zeros([size(L,1) 1]);zeros([size(D,1) 1])];
clear L D
%x=lsqr(H,y,1e-6,800 );
x=lsqr(H,y,1e-4,100,[],[],x0);

%fid=fopen('/RFrazinMatrices/xdyn_c2.14iml1e-5m100','wb')
fid=fopen('/RFrazinMatrices/xdyn_c2.14iml1e-5m1','wb');
fwrite(fid,x(:),'float32');
fclose(fid);

[x,rad,lat,lon,time] = dynamic_solution_movie('xdyn_c2.14iml1e-5m100','c2.14im','./RFrazinMatrices/',60,60,2.3,8);

 
x=reshape(x,60,60,120,14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% work on model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_fname='xmodel_c2'
%[x0,rad,lat,lon]=solution_movie(model_fname,directory,60,60,2.3,8.);
%solution2ascii(directory,model_fname,60,60,2.3,8.,6.);
[x0,rad,lat,lon] = readtom_sph(model_fname,directory,60,60,2.3,8.);
[A,y]=staticsparseread(directory,base,432000);
x0=x0(:);
y = A*x0(:);
% add noise
%yb = y +  .05*mean(y)*randn(size(y));
yb = y +  .02*mean(y)*randn(size(y));
% error functiom
norm_x0 = norm(x0,'fro');
err = @(x) norm(vec(x)-vec(x0))/norm_x0;

% unregularized solution
[x,flag,relres,iter,resvec,relvec] = lsqr(A,yb,[],200);
% relres=norm(y-A*x)/norm(y)=1.4 %
irad6 = 39;
snr(x,x0,irad6) % 0.57
solution_fname='xmodelrec_c2'; % without regularization

Rt0 = zeros(0,size(A,2)); % empty time regul matrix
Rs0 = zeros(0,size(A,2)); % empty spatial regul matrix
%build positive unregularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,y,Rs0,Rt0,0,0,xini,opts,1e-5);
solution_fname='xmodelrec_c2_nn'; % without regularization >0

lambda =1e-5;
[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
snr(x,x0,irad6) % 2.9
solution_fname='xmodelrec_c2_l1e-5';

%build positive regularized solution using TFOCS
xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,L,Rt0,lambda,0,xini,opts,1e-5);
snr(x,x0,irad6)  % 2.9
solution_fname='xmodelrec_c2_l1e-5_nn';

% gcv to compute solution%     
[m,n]=size(A);
u = randi(2,m,1);
u = 2*u-3; % -1,+1 random vector
score = @(q) gcvcalc([1e-5*q,0],A,yb,L,Rt0,u);
lambdamin=1e-5 ; lambdamax=1e2;
regparm_optim = fminbnd(score,lambdamin,lambdamax,...
    optimset('Display','iter'));
lambda = regparm_optim*1e-5;

lambda =2.5680e-7;% gcv solution (without noise)
[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
snr(x,x0,irad6) % 2.2
solution_fname='xmodelrec_c2_l2.6e-7'; %l=2.568e-7 gcv solution (without noise)

xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,L,Rt0,lambda,0,xini,opts,1e-5);
snr(x,x0,irad6) % 3
solution_fname='xmodelrec_c2_l2.6e-7_nn'

lambda =8.5217e-7;% gcv solution (with 5% noise)
[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
snr(x,x0,irad6) % 4.4 (2.5)
solution_fname='xmodelrec_c2_l8.6e-7'; %l=8.5217e-7 gcv solution (with 5% noise)

lambda =7.3047e-7;% gcv solution (77 images with 2% noise)
[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
snr(x,x0,irad6) % 3.56
solution_fname='xmodelrec_c2_77_l7.3e-7'; %l=7.3047e-7 gcv solution (77 images with 2% noise)

xini = x(:);
opts.tol=1e-4;
[x,out,opts] = solve_l2nonneg(A,yb,L,Rt0,lambda,0,xini,opts,1e-5);
snr(x,x0,irad6) % 4.85  (2.6257)  -> 14 im 
                % 3.57 -> 77 im
solution_fname='xmodelrec_c2_l8.6e-7_nn';
solution_fname='xmodelrec_c2_77_l7.3e-7nn';

% TV(theta,phi) regularized solution using TFOCS
clear opts;
opts.tol=1e-4;
opts.errFcn     = @(f,dual,primal) err(primal);
%epsilon = .05*mean(y)*sqrt(254040); % sigma*sqrt(n) ~ 30
epsilon = .02*mean(yb)*sqrt(1378196); % sigma*sqrt(n) ~ 30
nrad=60; ntheta=60; nphi=120;
%xini = x(:); % l=8.5217e-7 gcv solution (with 5% noise) nn
xini = A'*yb;
%xini = xunreg; % lsqr unregularized solution
An = 1e-5; % norm of A
[x,out,optsout] = solve_TVthetaphi(A,yb,epsilon,nrad,ntheta,nphi,xini,opts,An);
snr(x,x0,irad6) % 2.06 -> 14 im
                % 2.31 -> 77 im
solution_fname='xmodelrec_c2_TV_nn'
solution_fname='xmodelrec_c2_77_TV_nn'


solution_fname='xmodelrec_c2_TV_nn'

% write & display
fid=fopen([directory,solution_fname],'wb');
fwrite(fid,x(:),'float32');
fclose(fid);
[x,rad,lat,lon]=solution_movie(solution_fname,directory,60,60,2.3,8.);
 
solution_movie_mem(max(x,zeros(size(x))),rad,lon,lat);



% read solution
[x,rad,lat,lon] = readtom_sph(solution_fname,directory,60,60,2.3,8.);
x=x(:);

%check TV operator
TV=W(vec(x0),1);
solution_fname='xmodel_TV'
