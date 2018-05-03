function [x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,x0,opts,An)

    if nargin < 7, x0 = []; end
    if nargin < 8, opts = []; end
    if ~isfield( opts, 'restart' ), 
        opts.restart = 100; 
    end
    if nargin <9, An=linop_normest(A); end
        
    [s,n] = size(Rt);
    [t,n] = size(Rs);
    
    Areg = @(x,mode) Aregmult(x,mode,A/An,(lambda/An)*Rs,(mu/An)*Rt); 
    
    yfull = [y;zeros(s,1);zeros(t,1)];  
    
    [x,out,opts] = tfocs( smooth_quad, { Areg, -yfull }, proj_Rplus(), x0*An, opts );
    
    %scale back
    x = x/An;
end
