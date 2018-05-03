function x_cart = z_slice(x3D,n,sz)
%function im = z_slice(x3D,n,sz)
%
% this takes a slice at z = n and put it in a cartesian array
%    this is useful for displaying the slices
%
% x is a 3D array cylindrical array with dimensions
%   (nrad , nphi , nz)
% n  is the location of the z slice (image is taken at z = n)
% sz is the size of the sz by sz box in which the image
%    is to be displayed sz > 2*nrad and sz > 2*nphi
% x_cart is the output slice
%
%  see also ~rfrazin/ML/cart_to_cyl.m
%

dim = size(x3D); 
nrad = dim(1); nphi = dim(2); nz = dim(3);

if ( (sz < 2*nrad) | (sz < 2*nphi) )
  disp('SZ parameter is too small!!!')
  stop
end
if (n > nz)
  disp('bad number for n!!!')
  stop
end

cartsize = sz;

x_cyl2 = x3D(:,:,n);
x = zeros(sz,sz);
y = x;
x_cart  = x;

for i = 1:sz
  x(:,i) = i*ones(sz,1);
  y(i,:) = i*ones(1,sz);
end

x = x - mean(mean(x));
y = y - mean(mean(y));

r = sqrt(x.*x + y.*y);
r = .99999*r/max(max(x));
phi = atan2(y,x);
    ffi = find(phi < 0);
    phi(ffi) = phi(ffi) + 2*pi*ones(size(ffi));
    clear ffi x y

radind = 0;
phiind = 0;

for i = 1:sz
  for j = 1:sz
     radind = floor(   r(j,i)*nrad       ) + 1;
     phiind = floor( phi(j,i)*nphi/2./pi ) + 1;

     if (radind <= nrad)
        x_cart(j,i) = x_cyl2(radind,phiind );
     end 

   end
end