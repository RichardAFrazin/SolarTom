function solution2ascii(path,solution_fname,nrad,ntheta,Rmin,Rmax,Rmax_usefull);
%% convert to glnemo


endianstring = 'ieee-le';
[Ne,rad,lat,lon] = readtom_sph(solution_fname,path,nrad,ntheta,Rmin,Rmax,endianstring);
%[Ne,rad,lat,lon,time] = dynamic_readtom_sph(solution_fname,baseproj,path,nrad,ntheta,Rmin,Rmax,endianstring);
%ntime=length(time);
%Ne = Ne(:,:,:,1);

nnrad=find(rad > Rmax_usefull,1)-1;
rad=rad(1:nnrad); %throw beyond Rmax_usefull
Ne=Ne(1:nnrad,:,:);
Ne(Ne <= 0)=min(Ne(Ne > 0));

nspace=nnrad*ntheta^2*2;

dr=rad(2)-rad(1);
dtheta=abs(lat(2)-lat(1))*pi/180.;
dphi=abs(lon(2)-lon(1))*pi/180.;


%[lat,rad,lon] = meshgrid(lat,rad,lon);
[rad,lat,lon] = ndgrid(rad,lat,lon);
x = rad.*cosd(lat).*cosd(lon);
y = rad.*cosd(lat).*sind(lon);
z = rad.*sind(lat);

Volume = rad.^2.*cosd(lat)*dr*dtheta*dphi;
vsize = rad.^2*dr*dtheta*dphi;
size = 4*nthroot(3*vsize/(4*pi),3);

fid=fopen([path,solution_fname,'.txt'],'wb');

for i=1:nspace
    fprintf(fid,'%f %f %f %f %g %g\n',x(i),y(i),z(i),size(i),Ne(i),Volume(i)*Ne(i));    
end

fclose(fid);

% fid=fopen([path,solution_fname,'.txt'],'wb');
% for j=1:2*ntheta
%     for i=1:ntheta    
%         x = rad(1)*cosd(lat(i))*cosd(lon(j));
%         y = rad(1)*cosd(lat(i))*sind(lon(j));
%         z = rad(1)*sind(lat(i));
%         size = nthroot(rad(1)^2*dr*dtheta*dphi,3);
%         fprintf(fid,'%f %f %f %f %g\n',x,y,z,size,Ne(1,i,j,1));  
%     end
% end
% fclose(fid);
