function [d2r3d,d2theta3d,d2phi3d,h_laplac,r2_laplac,L] = spatialsparsebuildreg(nrad,ntheta,nphi,Rmin,Rmax)
% same as derivs_hollowsph.m
%
% D.VIBERT 20.12.2011

% index spherical coords
rad = linspace(Rmin + (Rmax-Rmin)/nrad/2, Rmax - (Rmax-Rmin)/nrad/2,nrad);
lat = linspace(-90 + 90/ntheta,90 - 90/ntheta,ntheta);
%lon = linspace(360/nphi/2, 360 - 360/nphi/2, nphi );

% truncated I matrices needed for 3d expansion via kronecker products
Ir     = sparse(1:nrad-2,2:nrad-1,ones(1,nrad-2),nrad-2,nrad);
Itheta = sparse(1:ntheta-2,2:ntheta-1,ones(1,ntheta-2),ntheta-2,ntheta);
Iphi   = speye(nphi);

% r derivatives 
dr = get_derivatives(nrad,1); % 1D 1st derivative 
dr = dr(2:nrad-1,:); % remove 1st line (to have same dim than 2nd derivative)
r = sparse(1:nrad-2,1:nrad-2,2*rad(2:nrad-1));
rdr = r*dr;
rdr3d = kron(Iphi,kron(Itheta,rdr)); % expand to 3D

d2r = get_derivatives(nrad,2);      % 1D 2nd derivative
r2 = sparse(1:nrad-2,1:nrad-2,rad(2:nrad-1).^2);% diag r^2
r2d2r = r2*d2r;
d2r3d   = kron(Iphi,kron(Itheta,d2r));  % expand to 3D
r2d2r3d = kron(Iphi,kron(Itheta,r2d2r));  % expand to 3D

% theta derivatives
dtheta = get_derivatives(ntheta,1); % 1D 1st derivative
dtheta = dtheta(2:ntheta-1,:); % remove 1st line (to have same dim than 2nd derivative)
sinocostheta = -sind(lat(2:ntheta-1))./cosd(lat(2:ntheta-1));
sinocostheta = sparse(1:ntheta-2,1:ntheta-2,sinocostheta);%diag -sin(theta)/cos(theta)
dtheta = sinocostheta*dtheta;
dtheta3d = kron(Iphi,kron(dtheta,Ir));  % expand to 3D

d2theta = get_derivatives(ntheta,2); % 1D 2nd derivative
d2theta3d = kron(Iphi,kron(d2theta,Ir)); % expand to 3D

% phi derivatives
d2phi = get_derivatives(nphi,2,1); % 1D 2nd derivative, periodic coordinate
oocos2lat = 1./(cosd(lat(2:ntheta-1)).^2);
oocos2lat = sparse(1:ntheta-2,2:ntheta-1,oocos2lat,ntheta-2,ntheta); % diag 1/cos(theta)^2
d2phi3d = kron(d2phi,kron(oocos2lat,Ir)); % expand to 3D

% true r^2 laplacian 
% r2laplac = r2 d2/dr2 + 2r d/dr + d2/d2theta - sin(theta)/cos(theta) d/dtheta
% +1/cos(theta)^2 d2/d2phi
%r2_laplac = r2d2r3d + rdr3d + d2theta3d + dtheta3d + d2phi3d ; 

% false laplacian
r2_laplac = r2d2r3d + d2theta3d  + d2phi3d;

% angular laplacian (remove constraint on d/dr)
% h_laplac = d2/d2theta + 1/cos(theta)^2 d2/d2phi
Ir     = speye(nrad);
d2theta3d = kron(Iphi,kron(d2theta,Ir)); % expand to 3D
d2phi3d   = kron(d2phi,kron(oocos2lat,Ir)); % expand to 3D
h_laplac = [d2theta3d ;d2phi3d]; 

% stacked derivatives 
Ir     = speye(nrad);
Itheta = speye(ntheta);
%oocos2lat = 1./(cosd(lat(1:ntheta)).^2);
%oocos2lat = sparse(1:ntheta,1:ntheta,oocos2lat); % diag 1/cos(theta)^2
%d2phi3d   = kron(d2phi,kron(oocos2lat,Ir)); % expand to 3D
d2phi3d   = kron(d2phi,kron(Itheta,Ir)); % expand to 3D
d2theta3d = kron(Iphi,kron(d2theta,Ir)); % expand to 3D
L = [d2theta3d ;d2phi3d]; 

return;