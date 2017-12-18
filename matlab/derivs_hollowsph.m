function [d2r,d2theta,d2phi,h_laplac] = derivs_hollowsph(nrad,ntheta,nphi,directory,fname_ext)
%function [d2r,d2theta,d2phi,h_laplac] = derivs_hollowsph(nrad,ntheta,nphi,directory,fname_ext)
%  this calculates the d2r,d2theta,d2phi derivs but returns the 
%     hlaplac = [d2theta;d2phi]  which is sort of a horizontal laplacian
%        only hlaplac is written to disk  in this version
%
%  Note that the final filenames will be will be *['hlaplac_',fname_ext]
%  
%
%  geometry is a hollow sphere grid with nz bins in z, nrad radial bins
%    and nphi angular bins.
%  each BLOCK is constant in Z
%  each COLUMN is constant in PHI
%   index =  (zind - 1)*nrows*ncols + (phiind - 1)*nrows + row;
%
%  it writes the matrices in ROW sparse format to disk
%    in the directory (inclue final slash)
%    and with the filename_extension
%
%
% note: in order for the row format sparse matrices to come out
%       correctly, the matrices must be built in r,theta,phi order.
%       i.e., r is the innermost, theta the second and phi the third   

if (nargin ~= 5)
  disp(['Wrong number of input arguments!']);
  return;
end


lplac = 0; % write 'horizontal laplacian'
all_three = 1;  % write deriv matrices for the 3 coords

nbins = nrad*ntheta*nphi;
row_d2r = zeros(1,3*nbins); col_d2r = row_d2r; val_d2r = row_d2r; 
row_d2phi = zeros(1,3*nbins); col_d2phi = row_d2phi; val_d2phi = row_d2phi;
row_d2theta = zeros(1,3*nbins); col_d2theta = row_d2theta; val_d2theta = row_d2theta;

%vectors of starting indecies for each row
nd2r = zeros(1,nbins+1); nd2phi = nd2r; nd2theta = nd2r;



%% calculate d2r

count = 0;
r_row_count = 0;
for k = 1:nphi
   for i = 1:ntheta
      for j = 2:nrad-1

            r_row_count = r_row_count + 1;
            
            count = count + 1;
            n = lindex3D(j-1,i,k,nrad,ntheta); %column
            row_d2r(count) = r_row_count;
            col_d2r(count) = n;
            val_d2r(count) = 1.0;

            count = count + 1;
            n = lindex3D(j,i,k,nrad,ntheta); %column
            row_d2r(count) = r_row_count;
            col_d2r(count) = n;
            val_d2r(count) = -2.0;

            count = count + 1;
            n = lindex3D(j+1,i,k,nrad,ntheta); %column
            row_d2r(count) = r_row_count;
            col_d2r(count) = n;
            val_d2r(count) = 1.0;
      
            nd2r(r_row_count+1) = count;
         
      end
   end
end

d2r  = sparse(row_d2r(1:count), col_d2r(1:count), val_d2r(1:count), r_row_count, nbins);
disp('done with d2r')

%calculate d2theta

count = 0;
t_row_count = 0;
for k = 1:nphi
   for i = 2:ntheta-1       
     for j = 1:nrad

        t_row_count = t_row_count + 1;

        count = count + 1;
        n = lindex3D(j,i,k,nrad,ntheta); %specify column
        row_d2theta(count) = t_row_count; %specify row
        col_d2theta(count) = n;
        val_d2theta(count) = -2.0;

        count = count + 1;
        n = lindex3D(j,i-1,k,nrad,ntheta); %column
        row_d2theta(count) = t_row_count;
        col_d2theta(count) = n;
        val_d2theta(count) = 1.0;

        count = count + 1;
        n = lindex3D(j,i+1,k,nrad,ntheta); %column
        row_d2theta(count) = t_row_count;
        col_d2theta(count) = n;
        val_d2theta(count) = 1.0;

        nd2theta(t_row_count+1) = count;

      end
   end
end

d2theta = sparse(row_d2theta(1:count),col_d2theta(1:count),val_d2theta(1:count),t_row_count,nbins);

val_hlaplac = val_d2theta(1:count);
col_hlaplac = col_d2theta(1:count);

disp('done with d2theta');

%% calculate d2phi

count = 0;
p_row_count = 0;
for k = 2:nphi-1
   for i = 1:ntheta
      for j = 1:nrad

            p_row_count = p_row_count + 1;

            count = count + 1;
            m = p_row_count; % specify row
            n = lindex3D(j,i,k-1,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;


            count = count + 1;
            n = lindex3D(j,i,k,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = -2.0;

            count = count + 1;
            n = lindex3D(j,i,k+1,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;      

            nd2phi(p_row_count+1) = count;
    
      end
   end
end
k = nphi;
   for i = 1:ntheta
      for j = 1:nrad

            p_row_count = p_row_count + 1;

            count = count + 1;
            m = p_row_count; % specify row
            n = lindex3D(j,i,k-1,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;


            count = count + 1;
            n = lindex3D(j,i,k,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = -2.0;

            count = count + 1;
            n = lindex3D(j,i,1,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;      

            nd2phi(p_row_count+1) = count;
    
      end
   end
k = 1;
   for i = 1:ntheta
      for j = 1:nrad

            p_row_count = p_row_count + 1;

            count = count + 1;
            m = p_row_count; % specify row
            n = lindex3D(j,i,nphi,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;


            count = count + 1;
            n = lindex3D(j,i,k,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = -2.0;

            count = count + 1;
            n = lindex3D(j,i,k+1,nrad,ntheta); % specify column
	    row_d2phi(count) = m;
            col_d2phi(count) = n;
            val_d2phi(count) = 1.0;      

            nd2phi(p_row_count+1) = count;
    
      end
   end


d2phi = sparse(row_d2phi(1:count), col_d2phi(1:count), val_d2phi(1:count), p_row_count,nbins);
disp('done with d2phi')

val_hlaplac = [val_hlaplac, val_d2phi(1:count)];
col_hlaplac = [col_hlaplac, col_d2phi(1:count)];
n_hlaplac = [nd2theta(1:t_row_count+1), nd2theta(t_row_count+1) + nd2phi(2:p_row_count+1)]; 

h_laplac = [d2theta;d2phi];
szh = size(h_laplac);

disp(['h_laplac has ',num2str(szh(1)), ...
      ' rows and ',num2str(szh(2)),' columns']);
 
disp(['output directory: ', directory])

if lplac 

   fname_hlaplac = ['hlaplac_',fname_ext];
   disp(['The filename extension is ', fname_hlaplac])

   y = zeros(t_row_count + p_row_count,1);

   fidn = fopen([directory,'n',fname_hlaplac],'wb');
   fidi = fopen([directory,'i',fname_hlaplac],'wb');
   fidv = fopen([directory,'v',fname_hlaplac],'wb');
   fidy = fopen([directory,'y',fname_hlaplac],'wb');
   fidd = fopen([directory,'delta_',fname_hlaplac],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_hlaplac]);
     disp([directory,'i',fname_hlaplac]);
     disp([directory,'v',fname_hlaplac]);
     disp([directory,'y',fname_hlaplac]);
	 disp([directory,'delta_',fname_hlaplac]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,n_hlaplac,'int32');
   fwrite(fidi,col_hlaplac-1,'int32');
   fwrite(fidv,val_hlaplac,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);

end %if laplac  



if all_three

   fname_d2r = ['d2r',fname_ext]

   y = zeros(r_row_count,1);
   
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


   fname_d2theta = ['d2theta',fname_ext]

   y = zeros(t_row_count,1);
   
   fidn = fopen([directory,'n',fname_d2theta],'wb');
   fidi = fopen([directory,'i',fname_d2theta],'wb');
   fidv = fopen([directory,'v',fname_d2theta],'wb');
   fidy = fopen([directory,'y',fname_d2theta],'wb');
   fidd = fopen([directory,'delta_',fname_d2theta],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_d2theta]);
     disp([directory,'i',fname_d2theta]);
     disp([directory,'v',fname_d2theta]);
     disp([directory,'y',fname_d2theta]);
     disp([directory,'delta_',fname_d2r]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,nd2theta,'int32');
   fwrite(fidi,col_d2theta-1,'int32');
   fwrite(fidv,val_d2theta,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);

   fname_d2phi = ['d2phi',fname_ext]

   y = zeros(p_row_count,1);

   fidn = fopen([directory,'n',fname_d2phi],'wb');
   fidi = fopen([directory,'i',fname_d2phi],'wb');
   fidv = fopen([directory,'v',fname_d2phi],'wb');
   fidy = fopen([directory,'y',fname_d2phi],'wb');
   fidd = fopen([directory,'delta_',fname_d2phi],'wb');

   if (fidn == -1 | fidi == -1 | fidv == -1 | ...
     fidy == -1 | fidd == -1 )
     disp('Names:');
     disp([directory,'n',fname_d2phi]);
     disp([directory,'i',fname_d2phi]);
     disp([directory,'v',fname_d2phi]);
     disp([directory,'y',fname_d2phi]);
     disp([directory,'delta_',fname_d2r]);
     error('Error opening one of the above files for writing');
   end
   
   fwrite(fidn,nd2phi,'int32');
   fwrite(fidi,col_d2phi-1,'int32');
   fwrite(fidv,val_d2phi,'float32');
   fwrite(fidy,y,'float32');
   fwrite(fidd,y,'float32');
 
   fclose(fidn);
   fclose(fidi);
   fclose(fidv);
   fclose(fidy);
   fclose(fidd);

end %all three

disp('done')
disp('')
disp('RUN row_to_col.c!!!!!')
disp('')
