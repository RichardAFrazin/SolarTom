Notes from email w/ Didier, 9/27/12
==========================

>[x,flag,res,iter] = solve_reg_lsqr(A,yb,L,Rt0,lambda,0,200);
>
>Est-ce que c'est la ligne pour résoudre le problème tomographique?

yes it is, for static tomo (since RT0 is empty and mu is zero)

the general use for dynamic tomo is
[x,flag,relres,iter,resvec] = solve_reg_lsqr(A,y,Rs,Rt,lambda,mu,200);
where,
 - Rs is the spatial regularization operator (applied to a 4D object)
 - Rt is the temporal regul op (on 4D)
The trick with Rt0 empty is you can use it for static tomo on 3D object

You can find some examples of how to use it for
 - static tomo: in static3Dsolve.m
 - dynamic tomo: in dynamic3Dsolve.m

You can also use the line
[x,out,opts] = solve_l2nonneg(A,y,Rs,Rt,lambda,mu,xini,opts,1e-5);
(see lines 105-110 in dynamic3Dsolve.m)
this is using TFOCS for positive constraint on L2 regularized problem
(no much longer to run than the non constrained version using lsmr algo)






