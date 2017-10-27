%this reads in the binary output file created by buildrow.c when
%  it is driven by geomtest.c

fid = fopen('geomtest.dat','rb');
aa = fread(fid,inf,'float64'); fclose(fid);
laa = length(aa);

%for 2 things in the array:
ind1 = [1:2:length(aa)]; ind2 = ind1 + 1;

q1 = aa(ind1); q2 = aa(ind2);