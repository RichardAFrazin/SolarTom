function [x,out,opts] = solve_TVthetaphi(A,y,epsilon,nrad,ntheta,nphi,x0,opts,An)

%{
    solve total-variation problem

    min_x ||x||_TVthetaphi
s.t.
    || A(x) - b || <= eps
    x>0

The solvers solve a regularized version


%}

% Before running this, please add the TFOCS base directory to your path



%% Call the TFOCS solver
opts.restart    = 1000;
%opts.maxIts     = 5000;
opts.maxIts     = 100;

clear tv;

W   = linop_TVthetaphi( [nrad,ntheta,nphi] );
normW      = linop_TVthetaphi( [nrad,ntheta,nphi], 'norm' );
normW2     = normW^2;
normA2     = An^2;
%mu = .005*norm( W(vec(x0),1) ,Inf);
mu =.005*norm( W(vec(x0),1) ,Inf)/normW*An;
opts.continuation = false;
continuationOptions.maxIts=4;

z0  = [];   % we don't have a good guess for the dual
tic;

%% next comes from solver_sBPDN_W
% with a positivity constraint added
proxScale   = sqrt( normW2 / normA2 );
prox        = { prox_l2( epsilon ), proj_linf(proxScale) };
W           = linop_compose( W, 1 / proxScale );

[ x, out, optsOut ] =   tfocs_SCD(  proj_Rplus(), { A, -y; W, 0 }, prox, mu, vec(x0), z0, opts, continuationOptions );

time_TFOCS = toc;

end

