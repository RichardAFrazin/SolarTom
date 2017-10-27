function [score,x]=cvcalc(regparm,Ad,yd,As,ys,Rs,Rt)
    global current_x;
    
    current_x=[];
    
    lambda = regparm(1);
    mu     = regparm(2);
    
    niter=50;
    %[x,flag] = solve_reg_lsqr(Ad,yd,Rs,Rt,lambda,mu,niter,current_x);
    [x,flag] = solve_reg_lsqr(Ad,yd,Rs,Rt,lambda,mu,niter);
    score = sum((As*x-ys).^2);
    %current_x = x;
    
%    disp(['cvcalc(',num2str(lambda(1)),',',num2str(lambda(2)),')=',num2str(score)]);
    return
    
