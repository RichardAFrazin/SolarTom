function index = lindex3D(row,col,block,nrows,ncols)
%function index = lindex3D(row,col,block,nrows,ncols)
%  this gives the 1D index of the (row,col,block) element
%    of a 3D matrix with nrows and ncols
%
%   see also vecdex3D.m

index = (block - 1)*nrows*ncols + (col - 1)*nrows + row;