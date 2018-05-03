function sphere_ascii(outfile, xname, r0, n_long, n_lat)
%function sphere_ascii(outfile, xname, r0, n_long, n_lat)

long = linspace(0, 360 - 360 / n_long, n_long);
lat = linspace(-90 + 180 / n_lat, 90 - 180 / n_lat , n_lat);

fid = fopen(outfile, 'w');
if (fid == -1)
  error(['Could not open file: ', outfile]);
end

[nrad, nphi, nz, rmin, rmax, grid_rmax, bindir, x_infile, cv_infile, ...
 imsize, build, type, outfile ] = get_build_opts;

[x_data, r_data, phi_data, z_data, nrad, nphi, nz] = readtom_cyl(xname);

if (0)
X = zeros(n_lat, n_long);
end

for i = 1:n_long
  longitude = long(i)*pi/180;
      
  for j = 1:n_lat
    latitude = lat(j)*pi/180;
    
    z = r0*sin(latitude);
    x = r0*cos(latitude)*cos(longitude);
    y = r0*cos(latitude)*sin(longitude);
    
    phi = atan2(y,x);
    if (phi < 0)
      phi = phi + 2 * pi;
    end
    
    [i_rad,i_phi,i_z] = cyl2ind(sqrt(x^2 + y^2), phi, z);
      
    a = [longitude*180/pi, latitude*180/pi, x_data(i_rad, i_phi, i_z)];
    fprintf(fid, '%.2f\t\t%.2f\t\t%g\n', a);
    
    if(0)
    X(j,i) = x_data(i_rad, i_phi, i_z);
    end
  end
end

fclose(fid);