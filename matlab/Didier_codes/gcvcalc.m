function [score,x,normRes,normRegs,normRegt]=gcvcalc(regparm,A,y,Rs,Rt,u)
    global current_x;
    
    current_x=[];
    
    lambda = regparm(1);
    mu     = regparm(2);
    
    niter=50;
    %[x,flag] = solve_reg_lsqr(Ad,yd,Rs,Rt,lambda,mu,niter,current_x);
    [x,flag] = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,niter);
    
    normRes  = norm(A*x-y);  
    normRegs = norm(Rs*x);
    normRegt = norm(Rt*x);
    
    % unbiased estimator of trace(I-A*(A'*A+lambda*Rs+mu*Rt)^-1*A')
    [m,n]=size(A);
%     u = randi(2,m,1);
%     u = 2*u-3; % -1,+1 random vector
    [xu,flag] = solve_reg_lsqr(A,u,Rs,Rt,lambda,mu,niter);
    denom = sum(u.*(u - A*xu));
    
    score = m*normRes^2/denom^2;
    %current_x = x;
    
%    disp(['cvcalc(',num2str(lambda(1)),',',num2str(lambda(2)),')=',num2str(score)]);
    return
    
