function [H,y,r] = dynamicsparsebuildproj(directory,base,objectsize)
%directory is where the files are, incl. final slash
%base is the suffix for the filename base 
%objectsize is the number of elements of the tomographic grid.  see
%   file ['info_',base]
%numimages is the number of images in the data (ie number of timesteps)
%
% D.Vibert & R.Frazin 7/9/2011

w = binfileread([directory,'w'],base,'float32'); %matrix entries
y = binfileread([directory,'y'],base,'float32'); %data vector
m = binfileread([directory,'m'],base,'int32');   %row start indicies
j = binfileread([directory,'j'],base,'int32');   %col indicies

r=binfileread([directory,'block_'],base,'int32'); %submatrices indicies
numimages = r(1);
r = r(2:numimages+2);

%to use the matlab sparse functions we need to generate the row 
%    indicies
j = j + 1; %matlab index convention
m = m + 1;
r = r + 1;
%rsteps = circshift(r,-1)-r;

k = zeros(size(w)); %row index vector

for l = 1:length(m)-1
    k(m(l):m(l+1)-1) = l;
end
%k(m(l+1):length(w)) = l+1;

% dynamic shape
% move the colums to the right place for each submatrices
for l=2:length(r)-1 % start with the 2nd submatrix
    nstart = m(r(l));
    nend   = m(r(l+1))-1;
    j(nstart:nend) = j(nstart:nend) + objectsize*(l-1);
end

H = sparse(k,j,w,length(m)-1,objectsize*numimages);

return;

