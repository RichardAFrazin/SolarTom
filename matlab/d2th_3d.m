function [d2th] = d2th_3d(m,n,o)
% this calculates the 
%    (1/r^2)*d2/d(theta)^2
%   matrix for an m X n X o object
% m = number of rows (y), n = number of cols (x
%

if (nargin ~= 3)
  error('delsq_sph: must have 3 input arguments')
end

mno = m*n*o;
mn = m*n;

[dx,dy,dz] = del_cen3d(m,n,o);
[d2x,d2y,d2z] = delsq3d(m,n,o);

   %the additional .5 is because this uses forward
   % differences 

x = zeros(m,n,o);
for i = 1:n
  x(:,i,:) = i*ones(m,1,o);
end
x = reshape(x,mno,1);
x = x - mean(x);

y = zeros(m,n,o);
for i = 1:m
  y(i,:,:) = i*ones(1,n,o);
end
y = reshape(y,mno,1);
y = y - mean(y);

z = zeros(m,n,o);
for i = 1:o
  z(:,:,i) = i*ones(m,n,1);
end
z = reshape(z,mno,1);
z = z - mean(z);

r = sqrt(x.*x + y.*y + z.*z) + 10*eps*ones(mno,1);

phi = sparse(atan2(y, ...
        x + 10*eps*ones(mno,1)));
theta = sparse(atan2(z, ...
        sqrt(x.*x + y.*y) + 10*eps*ones(mno,1)));

st = sin(theta); ct = cos(theta);
   clear theta
sp = sin(phi); cp =  cos(phi);
   clear phi

x = x./r; y = y./r; z = z./r;
   clear r
        

d2th =  2*diag(sparse(z.*z.*cp.*sp))*(dx*dy) + ...
       -2*diag(sparse(z.*ct.*cp   ))*(dx*dz) + ...
       -2*diag(sparse(z.*ct.*sp   ))*(dy*dz) + ...
          diag(sparse(z.*z.*cp.*cp))*(d2x)   + ...
          diag(sparse(z.*z.*sp.*sp))*(d2y)   + ...
          diag(sparse(ct.*ct      ))*(d2z) ;


          


