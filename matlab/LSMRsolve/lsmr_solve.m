fn = 'CR2099.171.20.90_DR0'; obsize = 20*90*180; lambda = 0.5;niter = 200;

[A,y]  = staticsparseread('/Users/frazin/tomography/bindata/',fn,obsize);
[Rs,ys] = staticsparseread('./','hlaplac_20_90',obsize);

nn = [1200,1600,2000];

lnn = length(nn);

for k = 1:lnn
    
  niter = nn(k);

  eval(['tic;[x',num2str(niter), ...
    ',flag',num2str(niter), ...
    ',iter',num2str(niter), ...
    ',relres',num2str(niter), ...
    '] = solve_reg_lsqr(A,y,Rs,lambda,niter); toc']);

end

