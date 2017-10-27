function [x,out,optsOut] = solve_TVthetaphi_TVtime(A,y,epsilon,alpha,beta,...
    nrad,ntheta,nphi,ntime,x0,opts,An,Rt)

%{
    solve total-variation problem

    min_x alpha ||x||_TVthetaphi + beta ||x||_TVtime
s.t.
    || A(x) - b || <= eps
  & x>0

The solvers solve a regularized version


%}

% Before running this, please add the TFOCS base directory to your path


if isempty(alpha), 
    alpha = 1; 
end
if isempty(beta), 
    beta = 1; 
end

%% Call the TFOCS solver
opts.restart    = 1000;
%opts.maxIts     = 5000;
opts.maxIts     = 150;

clear tv;

W   = linop_TVthetaphi4D( [nrad,ntheta,nphi,ntime] );
if isempty(Rt), 
    Wtime  = linop_dt4D([nrad*ntheta*nphi,ntime]);
else
    Wtime = Rt;
end

normW      = linop_TVthetaphi4D( [nrad,ntheta,nphi,ntime], 'norm' );

if isempty(Rt), 
    normWtime  = linop_dt4D( [nrad*ntheta*nphi,ntime], 'norm' );
else
    normWtime = linop_normest(Wtime);
end


opts.continuation = false;
continuationOptions.maxIts=4;

z0  = [];   % we don't have a good guess for the dual
tic;

%% next comes from solver_sBPDN_W
% with a positivity constraint added
proxScale1   = normW / An ;
proxScale2   = normWtime / An ;

prox        = { prox_l2( epsilon ), ...
                proj_linf( proxScale1 * alpha ),...
                proj_linf( proxScale2 * beta ) };
W1          = linop_compose( W, 1 / proxScale1 );
W2          = linop_compose( Wtime, 1 / proxScale2 );

mu =.005*norm( W1(vec(x0),1) ,Inf);

[ x, out, optsOut ] = tfocs_SCD(  proj_Rplus(), { A, -y; W1, 0; W2, 0 },...
                        prox, mu, vec(x0), z0, opts, continuationOptions );

time_TFOCS = toc;

end

