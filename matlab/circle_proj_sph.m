function [c_proj,y,z,lat,long] = circle_proj_sph(fname, r0, N, carlong, sun_ob1, pth)
%function [c_proj,y,z,lat,long] = circle_proj_sph(fname, r0, N,carlong, sun_ob1, [path])
%            
%   This projects a spherical shell of a tomographic reconstruction onto
%   the plane of the sky. 
% c_proj - N by N array of output image
% y,z - image coordinates (plane of sky)
% lat, long = latitude and longitude of image points
% fname - filename of object
% r0 - radius of spherical shell to be seen in projection (Rsun units)
% N - c_proj is N by N
% carlong - Carrington longitude of central meridian (in degrees)
% sun_ob1 = Sun --> spacecraft vector in J2000 cordinates (3x1 or 1x3)
%          (units don't matter)
%
%  note: this code does not include finite distance effects (assumes
%       parallel lines of sight)
%
%    see also circle_proj_cyl.m

%disp(['carlong = ',num2str(carlong),' degrees']);

%disp('circle_proj_sph: Mk4 parameters!');
%Rmax = 2.85; Rmin = 1.1; nrad = 35; ntheta = 60; nphi = 2*ntheta;
%disp([Rmax, Rmin, nrad, ntheta]);

%EIT/EUVI parameters
Rmax = 1.26; Rmin = 1.; 
%nrad = 13; ntheta = 60; 
nrad = 26; ntheta = 90;
nphi = 2*ntheta;

if (nargin ~= 6)
  pth = '/Users/frazin/tomography/bindata/';
end

[crap1,crap2] = size(sun_ob1);
if ( (crap2 == 3) && (crap1 == 1))
    sun_ob1 = sun_ob1';
elseif ( (crap2 == 3) && (crap1 == 1))
    ;
else
    error('sun_ob1 has wrong dimensions');
end
    


ALPHApo = 286.11 * pi/180.0;
DELTApo = 63.85  * pi/180.0;

spol1 = [cos(DELTApo)*cos(ALPHApo), cos(DELTApo)*sin(ALPHApo), sin(DELTApo)]';

x_data = readtom_sph(fname,pth,nrad,ntheta,Rmax,Rmin,'ieee-le');
 

c_proj = zeros(N,N);
lat = zeros(N,N);
long = zeros(N,N);

y = linspace(-r0, r0, N);
z = linspace(-r0, r0, N);
      
% compute R23 matrix
Rz = rotz( -atan2(sun_ob1(2), sun_ob1(1)));
sob = Rz*sun_ob1;
Ry = roty( -atan2(sob(3), sob(1)));
Rtmp = Ry*Rz;
spol2 = Rtmp*spol1;
Rx = rotx( atan2(spol2(2), spol2(3)));
R12 = Rx*Rtmp;

sun_ob2 = R12*sun_ob1;
spol2 = R12*spol1;

pang = atan2(spol2(1), spol2(3));
Ry = roty(pang);
Rz = rotz( carlong*pi/180 );
R23 = Rz*Ry;
      

for row = 1:N
  for col = 1:N
    z_r = z(row);
    y_c = y(col);
    
    if (sqrt(y_c^2 + z_r^2) <= r0)
      x = sqrt(r0^2 - y_c^2 - z_r^2);
      
      r_2 = [x; y_c; z_r];
      r_3 = R23 * r_2;
      
      rad = sqrt(r_3(1)^2 + r_3(2)^2 + r_3(3)^2);
      theta = atan( r_3(3) / sqrt(r_3(1)^2 + r_3(2)^2) ) + pi/2;
      phi = atan2(r_3(2),r_3(1));
      if (phi < 0)
        phi = phi + 2*pi;
      end
      long(row,col) = phi*180/pi;
      lat(row,col) = (theta - pi/2)*180/pi;
    
    
      r_i = floor( nrad*(rad - Rmin)/(Rmax - Rmin) ) + 1 ;
      theta_i = floor( theta*ntheta/pi ) + 1;
      phi_i = floor( nphi*phi/2/pi) + 1;
      
      
      if (r_i < 1 | r_i > nrad)
       disp([row, col, z_r, y_c, rad, r_i]);
       error('r ssertion failed');
      end
      if (theta_i < 1 | theta_i > ntheta)
       disp([row, col, z_r, y_c, theta, theta_i]);
       error('theta assertion failed');
      end
      if (phi_i < 1 | phi_i > nphi)
       disp([row, col, z_r, y_c, phi, phi_i]);
       error('phi assertion failed');
      end
      
      c_proj(row,col) = x_data(r_i, theta_i, phi_i);
    else
      cproj(row,col) = 1.0;
    end
  end
end