%this prints out the tomographic reconstructions for Chip

MK4=0;
C2_90 = 0;
c2_60_60 = 0;
c2_40_60 = 1;
EUVI_26_90 = 0;

%pth = '/Users/frazin/tomography/bindata/';
%fn = 'x_c2MARS.pb.2008.1204.1217.rmax8.0nr60nt60l3e-5';
pth = '/Users/frazin/Desktop/';
fn = 'x_c2_40_60_b3_2005.0430.0513.nocme.mb_l3e-5';


if MK4
  Rmin = 1.1; Rmax = 1.85; nrad = 40; ntheta = 60; endian = 'ieee-le';
end
if C2_90
  Rmin = 2.3; Rmax = 6.1; nrad = 38; ntheta = 90; endian = 'ieee-le';
end
if c2_40_60
  Rmin = 2.3; Rmax = 6.1; nrad = 40; ntheta = 60; endian = 'ieee-le';
end
if EUVI_26_90
  Rmin = 1.0; Rmax = 1.26; nrad = 26; ntheta = 90; endian = 'ieee-le';
end
if c2_60_60   
  Rmin = 2.3; Rmax = 8.0;  nrad = 60; ntheta = 60; endian = 'ieee-le';  
end
  
split_in_two = 0;  %see below and edit
cartesian = 1;
spherical = 0;  




[dat,rad,lat,lon] = readtom_sph(fn,pth,nrad,ntheta,Rmin,Rmax,endian);

lat = lat*pi/180; lon = lon*pi/180;

%for formats see 'doc fprintf'

%disp(size(dat));
%disp(size(rad));
%disp(size(lat));
%disp(size(lon));
%return;

if split_in_two 

  fid = fopen('/Users/frazin/Desktop/crap_pt1','w');
  for i = 1:nrad/2
   for j = 1:ntheta
    for k = 1:2*ntheta
      x = rad(i)*cos(lat(j))*cos(lon(k));
      y = rad(i)*cos(lat(j))*sin(lon(k));
      z = rad(i)*sin(lat(j));
      fprintf(fid,'%7.3f %7.3f %7.3f %6.2g\n',x,y,z,dat(i,j,k));
      %fprintf(fid,'%3.1f %4.1f %4.1f %6.2g\n',rad(i),lat(j),lon(k),dat(i,j,k));
    end
   end
  end
  fclose(fid);

  fid = fopen('/Users/frazin/Desktop/crap_pt2','w');
  for i = (nrad/2+1):nrad
   for j = 1:ntheta
    for k = 1:2*ntheta
      x = rad(i)*cos(lat(j))*cos(lon(k));
      y = rad(i)*cos(lat(j))*sin(lon(k));
      z = rad(i)*sin(lat(j));
      fprintf(fid,'%7.3f %7.3f %7.3f %6.2g\n',x,y,z,dat(i,j,k));
      %fprintf(fid,'%3.1f %4.1f %4.1f %6.2g\n',rad(i),lat(j),lon(k),dat(i,j,k));
    end
   end
  end
  fclose(fid);
  
else
  fid = fopen(['/Users/frazin/Desktop/',fn,'.txt'],'w');
  if c2_60_60
      iimax = 40;
  else
      iimax = nrad
  end
  
  if (cartesian == 1)
    for i = 1:iimax
     for j = 1:ntheta
      for k = 1:2*ntheta
        x = rad(i)*cos(lat(j))*cos(lon(k));
        y = rad(i)*cos(lat(j))*sin(lon(k));
        z = rad(i)*sin(lat(j));
        fprintf(fid,'%+7.6f %+7.6f %+7.6f %+6.3g\n',x,y,z,dat(i,j,k));
       %fprintf(fid,'%3.1f %4.1f %4.1f %6.2g\n',rad(i),lat(j),lon(k),dat(i,j,k));
      end
     end
    end
  elseif (spherical == 1)
    for i = 1:iimax
      for j = 1:ntheta
         for k = 1:2*ntheta
              fprintf(fid,'%+7.6f %+7.6f %+7.6f %+6.3g\n', ...
                   rad(i),lat(j),lon(k),dat(i,j,k));
         end
      end
    end
  end


    fclose(fid);

end
