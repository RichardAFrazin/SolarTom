function [c_proj,y,z,lat,long] = circle_proj_cyl(fname, r0, N, carlong, sun_ob1, path)
%function [c_proj,y,z,lat,long] = circle_proj_cyl(fname, r0, N,
%                                              carlong, sun_ob1, [path])
%   see also circle_proj_sph.m


ALPHApo = 286.11 * pi/180.0;
DELTApo = 63.85  * pi/180.0;

spol1 = [cos(DELTApo)*cos(ALPHApo), cos(DELTApo)*sin(ALPHApo), sin(DELTApo)]';

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

if (nargin == 6)
  [x_data, r_data, phi_data, z_data, nrad, nphi, nz] = ...
      readtom_cyl(fname, path);
else
  [x_data, r_data, phi_data, z_data, nrad, nphi, nz] = ...
      readtom_cyl(fname);
end

c_proj = zeros(N);
lat = zeros(N);
long = zeros(N);

z_min = -grid_rmax;
z_max = grid_rmax;

y = linspace(-r0, r0, N);
z = linspace(-r0, r0, N);

deltagrid = 2 * grid_rmax / nz;
      
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
      
      rad = sqrt(r_3(1)^2 + r_3(2)^2);
      phi = atan2(r_3(2), r_3(1));
      
      if (phi < 0)
        phi = phi + 2*pi;
      end
      
      long(row, col) = phi;
      lat(row, col) = atan2(r_3(3),sqrt(r_3(1)^2+r_3(2)^2));
      
      z_3 = r_3(3);
      
      z_i = floor((z_3 + grid_rmax) / deltagrid) + 1;
      
      r_i = floor((rad / grid_rmax) * nrad) + 1;
      phi_i = floor((phi / 2 / pi) * nphi) + 1;
      
      if (z_i < 1 | z_i > nz | ...
          r_i < 1 | r_i > nrad | ...
          phi_i < 1 | phi_i > nphi)
        error('Assertion failed');
      end
      
      c_proj(row,col) = x_data(r_i, phi_i, z_i);
    end
  end
end