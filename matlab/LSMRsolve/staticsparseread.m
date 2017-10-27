function [H,y] = staticsparseread(directory,base,objectsize)
%directory is where the files are, incl. final slash
%base is the suffix for the filename base 
%objectsize is the number of elements of the tomographic grid.  see
%   file ['info_',base]
%
% D.Vibert & R.Frazin 7/9/2011

w = binfileread([directory,'w'],base,'float32'); %matrix entries
y = binfileread([directory,'y'],base,'float32'); %data vector
m = binfileread([directory,'m'],base,'int32');   %row start indicies
j = binfileread([directory,'j'],base,'int32');   %col indicies

%to use the matlab sparse functions we need to generate the row 
%    indicies

j = j + 1; %matlab index convention
m = m + 1;

k = zeros(size(w)); %row index vector

for l = 1:length(m)-1
    k(m(l):m(l+1)-1) = l;
end
%k(m(l+1):length(w)) = l+1;


H = sparse(k,j,w,length(m)-1,objectsize);

return;

