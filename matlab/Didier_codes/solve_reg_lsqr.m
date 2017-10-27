function [x,flag,relres,iter] = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,niter)
%%
% A projection matrix (or function handle)
% y right hand side vector 
% Rs spatial regul matrix 
% Rt temporal regul matrix
% lambda spatial regul parameter
% mu time regul parameter
%
%% D.Vibert 16/12/2011 

    Areg = @(x,mode) Aregmult(x,mode,A,lambda*Rs,mu*Rt); 
    
    if nargin < 8 , x0=[];
        
    [s,n] = size(Rt);
    [t,n] = size(Rs);
    
    yfull = [y;zeros(s,1);zeros(t,1)];  
%    [x,flag,relres,iter,resvec,relvec] = lsqr(Areg,yfull,[],niter);

%   better algo from standford group (minimize faster the normal equation
%   residual)
    [x,flag,iter,relres] = lsmr(Areg,yfull,[],[],[],[],niter,[]);
    
return
end


