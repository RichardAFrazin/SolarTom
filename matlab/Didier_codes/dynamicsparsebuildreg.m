function [L,D] = dynamicsparsebuildreg(directory,baseproj,basereg)
%directory is where the files are, incl. final slash
%baseproj is the suffix for the filename base of the projections
%basereg                             of the spatial regularization matrix
%objectsize is the number of elements of the tomographic grid.  see
%   file ['info_',base]
%numimages is the number of images in the data (ie number of timesteps)
%
% D.Vibert & R.Frazin 9/9/2011 

w = binfileread([directory,'w'],basereg,'float32'); %matrix entries
%y = binfileread([directory,'y'],basereg,'float32'); %data vector
m = binfileread([directory,'m'],basereg,'int32');   %row start indicies
j = binfileread([directory,'j'],basereg,'int32');   %col indicies

% read dates in mjd
mjd = binfileread([directory,'date_'],baseproj,'float64');

% check read was ok
if ((length(w) == 1) | (m(1) == -1) | (j(1) == -1) | (mjd(1) == -1) )
    disp(['Problem reading files. w=',num2str(w), ...
        ' m=',num2str(m),' j=',num2str(j),' mjd=',num2str(mjd) ]);
end

numimages = length(mjd);
deltat = circshift(mjd,-1) - mjd;
deltat = deltat(1:numimages-1);

%to use the matlab sparse functions we need to generate the row 
%    indicies
j = j + 1; %matlab index convention
m = m + 1;
static_objectsize = max(j);

k = zeros(size(w)); %row index vector

for l = 1:length(m)-1
    k(m(l):m(l+1)-1) = l;
end

% one static block
LL = sparse(k,j,w,length(m)-1,static_objectsize);

% build full dynamic spatial regul matrix
L = kron(speye(numimages),LL);
clear LL;

% build one temporal block
DT = diag(ones(size(deltat))./deltat,1);
DT = DT(1:numimages-1,:);
DT = -circshift(DT,[0,-1]) + DT;
DT = sparse(DT);

% build full dynamic temporal regul matrix
D = kron(DT,speye(static_objectsize));
clear DT;

return;

