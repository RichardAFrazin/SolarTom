function [nB,iB,vB] = sparfix(B);
%function [nB,iB,vB] = sparfix(B);
%B(nf,np) = original matrix
%nB(1,np) = the number of nz elements of each column of B
%iB(nzmax,np) = row indecies of those elements
%vB(nzmax,np) = values of those elements

sB = size(B); 
nr = sB(1); nc = sB(2); 
clear sB;

for ci = 1:nc
  nB(ci) = nnz(B(:,ci));
end 

nzmax = max(nB);
iB =   zeros(nzmax,nc);
vB =   zeros(nzmax,nc);

for ci = 1:nc
  ind =  find(B(:,ci));
  iB(1:nB(ci),ci) = ind;
  vB(1:nB(ci),ci) = B(ind,ci);
end

%use C-style row indecies 
iB = iB - 1;