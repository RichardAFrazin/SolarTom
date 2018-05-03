function data_cart = pol2cart_data(nphi,r,data_pol,n,minval)

nrad = length(r);

if (nargin == 5)
  data_cart = minval*ones(n);
else
  data_cart = zeros(n);
end

rmax = r(end);
rmin = r(1);

x = zeros(n);
y = x;

index = (1:n) - mean(1:n);

x = ones(n,1) * index;
y = flipud(index') * ones(1,n);

r = sqrt(x.*x + y.*y);
r = r / max(max(r));
r = r * inv(r(1,floor(n/2))) * rmax;
r = r - rmin;
r = floor(r / r(1,floor(n/2)+1) * nrad);

phi = atan2(y,x) - pi/2;
I = find(phi < 0);
phi(I) = phi(I) + 2*pi*ones(size(I));

for i = 1:n
  for j = 1:n
    radind = r(j,i);
    phiind = floor(phi(j,i)*nphi/2/pi) + 1;
    
    if (radind <= nrad & radind >= 1)
      data_cart(j,i) = data_pol(radind,phiind);
    end 
  end
end