%this prints out a spherical shell cut of a tomographic reconstruction for Chip

MK4=0;
C2_60 = 1;


pth = '/Users/frazin/tomography/bindata/';
fn = 'x_c2_40_60_b3_2005.0514.0527_l3e-5';

if MK4
  Rmin = 1.1; Rmax = 1.85; nrad = 40; ntheta = 60; endian = 'ieee-le';
end
if C2_60
  Rmin = 2.3; Rmax = 6.1; nrad = 40; ntheta = 60; endian = 'ieee-le';
end

[dat,rad,lat,lon] = readtom_sph(fn,pth,nrad,ntheta,Rmin,Rmax,endian);

%for formats see 'doc fprintf'
%disp(size(dat));
%disp(size(rad));
%disp(size(lat));
%disp(size(lon));
%return;

cartesian = 0;
spherical = 1;  

i = 8; %3.022 Rs for the C2 40x60x180

  fid = fopen('/Users/frazin/Desktop/cr2029.dat_r3.02_sph.txt','w');
  if (cartesian == 1)

     lat = lat*pi/180; lon = lon*pi/180;
    
     for j = 1:ntheta
      for k = 1:2*ntheta
        x = rad(i)*cos(lat(j))*cos(lon(k));
        y = rad(i)*cos(lat(j))*sin(lon(k));
        z = rad(i)*sin(lat(j));
        fprintf(fid,'%+7.6f %+7.6f %+7.6f %+6.3g\n',x,y,z,dat(i,j,k));
       %fprintf(fid,'%3.1f %4.1f %4.1f %6.2g\n',rad(i),lat(j),lon(k),dat(i,j,k));
      end
     end
    
  elseif (spherical == 1)
    
      for j = 1:ntheta
         for k = 1:2*ntheta
%              fprintf(fid,'%+7.6f %+7.6f %+7.6f %+6.3g\n', ...
%                   rad(i),lat(j),lon(k),dat(i,j,k));
              fprintf(fid,'%+7.6f %+7.6f %+6.3g\n', ...
                   lat(j),lon(k),dat(i,j,k));
         end
      end
   
  end


    fclose(fid);
