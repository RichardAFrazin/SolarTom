function [slice,r,z] = phi_slice(fname, phi_p, path, nrad_p)
% function [slice,r,z] = phi_slice(fname, phi_p, path, nrad_p)

if (nargin >=3)
  [x, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname,path);
else
  [x, r, phi, z, nrad, nphi, nz] = readtom_cyl(fname);
end

phi = phi_p;

if (nargin == 4)
  max_r = max(r);
  r = linspace(-max_r, max_r, nrad_p);
  nrad = nrad_p
else
  max_r = max(r);
  nrad = 2*nrad;
  r = linspace(-max_r, max_r, nrad);
end

slice = zeros(nz, nrad);

for i=1:nz
  for j=1:nrad
    if (r(j) < 0)
      %-r(j)
      %mod(phi+180,360)
      %z(i)
      [i_rad,i_phi,i_z] = cyl2ind(-r(j), mod(phi+180,360)*pi/180, z(i));
    else
      [i_rad,i_phi,i_z] = cyl2ind(r(j), phi*pi/180, z(i));
    end
    %i_rad
    %i_phi
    x(i_rad, i_phi, i_z);
    slice(i,j) = x(i_rad, i_phi, i_z);
  end
end