function [x,flag,relres,iter] = solve_reg_lsqr(A,y,Rs,lambda,niter)
%%
% A projection matrix (or function handle)
% y right hand side vector 
% Rs spatial regul matrix 
% lambda spatial regul parameter
%
%% D.Vibert 16/12/2011 - modified R. Frazin 10/16/12

    Areg = @(x,mode) Aregmult(x,mode,A,lambda*Rs); 
    
    [s,n] = size(Rs);

    yfull = [y;zeros(s,1)];  


%default tola=tolb = 1.e-6
tola = 1.e-10;
tolb = tola;
    
%lsmr solver
 [x,flag,iter,relres] = lsmr(Areg,yfull,[],tola,tolb,[],niter,[]);
    
return
end


