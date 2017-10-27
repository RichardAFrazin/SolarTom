function [xvec] = boundary_draw(ti,tf,pi,pf,nr,ntheta);
%function [xvec] = boundary_draw(ti,tf,pi,pf,nr,ntheta);
%
%This function makes a vector reprenting the voxels of a regular
%  spherical grid with a total of nc3 = nr*ntheta*(2*ntheta) components,
%  in which all components are zero except for the box (in the first
%  radial bin) OUTSIDE of (theta,phi) (ti,tf,pi,pf).  In those, the 
%  value of the vector will be -1.
%
%Note that this code assumes the input indices start at 1, not 0.  If
%  they start at 0 (e.g., C or IDL indices) then the input indices
%  should have 1 added to them.

nc3 = nr*ntheta*(2*ntheta);

left_bdy = zeros(tf-ti+2,3); right_bdy = left_bdy;
top_bdy  = zeros(pf-pi+2,3); bot_bdy =  top_bdy;
xvec = zeros(nc3,1);

rad = 1;

%left,right boundary pixels:
len = tf - ti + 3;
for k = 0:len-1
    m = k + 1;
    left_bdy(m,:) =   [rad, ti - 1 + k,pi - 1];
    right_bdy(m,:) =  [rad, ti - 1 + k,pf + 1];
    xvec(lindex3D(left_bdy(m,1),left_bdy(m,2),left_bdy(m,3),nr,ntheta)) = -1.;
    xvec(lindex3D(right_bdy(m,1),right_bdy(m,2),right_bdy(m,3),nr,ntheta)) = -1.;
end

len = pf - pi + 3;
for k = 0:len-1
    m = k + 1;
    top_bdy(m,:) = [rad, tf + 1, pi - 1 + k];
    bot_bdy(m,:) = [rad, ti - 1, pi - 1 + k];
    xvec(lindex3D(top_bdy(m,1),top_bdy(m,2),top_bdy(m,3),nr,ntheta)) = -1.;
    xvec(lindex3D(bot_bdy(m,1),bot_bdy(m,2),bot_bdy(m,3),nr,ntheta)) = -1.;
end



%left_bdy
%right_bdy
%top_bdy
%bot_bdy

