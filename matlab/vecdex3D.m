function [row,col,block] = vecdex3D(lindex,nrows,ncols)
%function [row,col,block] = vecdex3D(lindex,nrows,ncols)
%   this returns the 3D index of the (1D) lindex_th element
%     of a 3D matrix with nrows and ncols
%
%     see also lindex3D.m

   n = lindex - 1;
   block = floor(n/(nrows*ncols));
   m = n - block*nrows*ncols;
   col = floor(m/nrows);
   row = m - nrows*col;
   block = block + 1;
   row = row + 1;
   col = col + 1;

