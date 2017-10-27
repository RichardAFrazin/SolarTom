function [d2r,dphi,d2z] = derivs_cyl_occ(nrad,nphi,nz,Rmin,grid_Rmax,directory,fname_ext)
%function [d2r,dphi,d2z] = derivs_cyl_occ(nrad,nphi,nz,Rmin,grid_Rmax,directory,fname_ext)
%  this calculates the d2r,dphi,d2z derivs
%  on a cylindrical grid with nz bins in z, nrad radial bins
%    and nphi angular bins.
%  each BLOCK is constant in Z
%  each COLUMN is constant in PHI
%   index =  (zind - 1)*nrows*ncols + (phiind - 1)*nrows + row;
%
%  Rmin is the occulted radius
%  grid_Rmax is the radius of the computation box
%  This puts ones for the central region (r <= Rmin) in the d2z matrix
%    and only regularizes outside of the central region in the d2r matrix
%    the dphi matrix doesnt change 
%
%  it writes the matrices in ROW sparse format to disk
%    in the directory (inclue final slash)
%    and with the filename_extension
%
%
% note: in order for the row format sparse matrices to come out
%       correctly, the matrices must be built in r,phi,z order.
%       i.e., r is the innermost, phi the second and z the third   

if (nargin ~= 7)
  disp(['Wrong number of input arguments!']);
  return;
end

nbins = nrad*nphi*nz;
row_d2r = zeros(1,3*nbins); col_d2r = row_d2r; val_d2r = row_d2r; 
row_d2z = zeros(1,3*nbins); col_d2z = row_d2z; val_d2z = row_d2z;
row_dphi = zeros(1,2*nbins); col_dphi = row_dphi; val_dphi = row_dphi;

%vectors of starting indecies for each row
nd2r = zeros(1,nbins+1); nd2z = nd2r; ndphi = nd2r;

rho = zeros(nrad,nphi,nz);
r = linspace(grid_Rmax/nrad,grid_Rmax,nrad);
z = linspace(-grid_Rmax,grid_Rmax,nz);
for k = 1:nz
   for i = 1:nphi
      for j = 1:nrad
         rho(j,i,k) = sqrt(r(j)^2 + z(k)^2);
      end
   end
end


%% calculate d2z

count = 0;
row_count = 0;
for k = 1:nz
   for i = 1:nphi
      for j = 1:nrad
       if ((k > 1) & (k < nz))

         rho1 = rho(j,i,k-1); rho2 = rho(j,i,k); rho3 = rho(j,i,k+1);
	 if ((rho1 > Rmin) & (rho3 > Rmin))

            count = count + 1;
            m = lindex3D(j,i,k,nrad,nphi); % specify row
            n = lindex3D(j,i,k-1,nrad,nphi); % specify column
	    row_d2z(count) = m;
            col_d2z(count) = n;
            val_d2z(count) = 1.0;


            count = count + 1;
            n = lindex3D(j,i,k,nrad,nphi); % specify column
	    row_d2z(count) = m;
            col_d2z(count) = n;
            val_d2z(count) = -2.0;

            count = count + 1;
            n = lindex3D(j,i,k+1,nrad,nphi); % specify column
	    row_d2z(count) = m;
            col_d2z(count) = n;
            val_d2z(count) = 1.0;      

         elseif (rho2 <= Rmin)

            count = count + 1;
            n = lindex3D(j,i,k,nrad,nphi); % specify column
	    row_d2z(count) = m;
            col_d2z(count) = n;
            val_d2z(count) = 2.0;

         end   
       end   

       row_count = row_count + 1;
       nd2z(row_count+1) = count;
    
      end
   end
end

d2z = sparse(row_d2z(1:count), col_d2z(1:count), val_d2z(1:count), nbins,nbins);
disp('done with d2z')

%% calculate d2r

count = 0;
row_count = 0;
for k = 1:nz
   for i = 1:nphi

      rho1 = rho(1,i,k);

      if ( rho1 > Rmin )
         count = count + 1;
         m = lindex3D(1,i,k,nrad,nphi); %specify row
         n = lindex3D(1,i,k,nrad,nphi); %column
         row_d2r(count) = m;
         col_d2r(count) = n;
         val_d2r(count) = -1.0;

         count = count + 1;
         n = lindex3D(2,i,k,nrad,nphi); %column
         row_d2r(count) = m;
         col_d2r(count) = n;
         val_d2r(count) = 1.0;
      end

      row_count = row_count + 1;
      nd2r(row_count+1) = count;
      

      for j = 2:nrad
       if (j < nrad)

         rho1 = rho(j-1,i,k);

         if (rho1 > Rmin)
            count = count + 1;
            m = lindex3D(j,i,k,nrad,nphi); %specify row
            n = lindex3D(j-1,i,k,nrad,nphi); %column
            row_d2r(count) = m;
            col_d2r(count) = n;
            val_d2r(count) = 1.0;

            count = count + 1;
            n = lindex3D(j,i,k,nrad,nphi); %column
            row_d2r(count) = m;
            col_d2r(count) = n;
            val_d2r(count) = -2.0;

            count = count + 1;
            n = lindex3D(j+1,i,k,nrad,nphi); %column
            row_d2r(count) = m;
            col_d2r(count) = n;
            val_d2r(count) = 1.0;
         end
       end
      
       row_count = row_count + 1;
       nd2r(row_count+1) = count;
         
      end
   end
end

d2r  = sparse(row_d2r(1:count), col_d2r(1:count), val_d2r(1:count), nbins, nbins);
disp('done with d2r')


count = 0;
row_count = 0;
for k = 1:nz

      for j = 1:nrad

      count = count + 1;
      m = lindex3D(j,1,k,nrad,nphi); %specify row
      n = lindex3D(j,1,k,nrad,nphi); %column
      row_dphi(count) = m;
      col_dphi(count) = n;
      val_dphi(count) = -1.0;

      count = count + 1;
      n = lindex3D(j,nphi,k,nrad,nphi); %column
      row_dphi(count) = m;
      col_dphi(count) = n;
      val_dphi(count) = 1.0;

      row_count = row_count + 1;
      ndphi(row_count+1) = count;

      end

   for i = 2:nphi
       
      for j = 1:nrad

      count = count + 1;
      m = lindex3D(j,i,k,nrad,nphi); %specify row
      n = lindex3D(j,i,k,nrad,nphi); %column
      row_dphi(count) = m;
      col_dphi(count) = n;
      val_dphi(count) = 1.0;

      count = count + 1;
      n = lindex3D(j,i-1,k,nrad,nphi); %column
      row_dphi(count) = m;
      col_dphi(count) = n;
      val_dphi(count) = -1.0;

      row_count = row_count + 1;
      ndphi(row_count+1) = count;

      end
   end
end

dphi = sparse(row_dphi(1:count),col_dphi(1:count),val_dphi(1:count),nbins,nbins);

disp('done with dphi')


   
   format compact

   y = zeros(nbins,1);

   fname_d2r = ['d2r',fname_ext]
   
   fidn = fopen([directory,'n',fname_d2r],'wb');
   fidi = fopen([directory,'i',fname_d2r],'wb');
   fidv = fopen([directory,'v',fname_d2r],'wb');
   fidy = fopen([directory,'y',fname_d2r],'wb');
   fidd = fopen([directory,'delta_',fname_d2r],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_d2r]);
     disp([directory,'i',fname_d2r]);
     disp([directory,'v',fname_d2r]);
     disp([directory,'y',fname_d2r]);
     disp([directory,'delta_',fname_d2r]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,nd2r,'int32');
   fwrite(fidi,col_d2r-1,'int32');
   fwrite(fidv,val_d2r,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);

   fname_d2z = ['d2z',fname_ext]
   
   fidn = fopen([directory,'n',fname_d2z],'wb');
   fidi = fopen([directory,'i',fname_d2z],'wb');
   fidv = fopen([directory,'v',fname_d2z],'wb');
   fidy = fopen([directory,'y',fname_d2z],'wb');
   fidd = fopen([directory,'delta_',fname_d2z],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_d2z]);
     disp([directory,'i',fname_d2z]);
     disp([directory,'v',fname_d2z]);
     disp([directory,'y',fname_d2z]);
     disp([directory,'delta_',fname_d2r]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,nd2z,'int32');
   fwrite(fidi,col_d2z-1,'int32');
   fwrite(fidv,val_d2z,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);

   fname_dphi = ['dphi',fname_ext]
   
   fidn = fopen([directory,'n',fname_dphi],'wb');
   fidi = fopen([directory,'i',fname_dphi],'wb');
   fidv = fopen([directory,'v',fname_dphi],'wb');
   fidy = fopen([directory,'y',fname_dphi],'wb');
   fidd = fopen([directory,'delta_',fname_dphi],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_dphi]);
     disp([directory,'i',fname_dphi]);
     disp([directory,'v',fname_dphi]);
     disp([directory,'y',fname_dphi]);
     disp([directory,'delta_',fname_d2r]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,ndphi,'int32');
   fwrite(fidi,col_dphi-1,'int32');
   fwrite(fidv,val_dphi,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);



disp('done')
disp('')
disp('RUN row_to_col.c!!!!!')
disp('')
