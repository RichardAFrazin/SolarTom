function [a,ah] = read_mars(filename,ppath,r_block)
%function read_mars(filename,r_block)
%   This reads the Marseilles LASCO pB images.
%a - 512x512 output image with blocking
%ah - fits header structure
%filename  - duh
%ppath - directory where filename resides.  include final slash
%r_block - optional argument for the blocking radius in units of
%   Rsun where Rsun is taken to be the instaneously correct value. Default
%   value is 2.3

blocked_pix_value = .1;

if ((nargin < 2) || (nargin > 3))
    disp('read_mars.m: Wrong number of input arguments!');
    return;
end
if nargin == 2
    r_block = 2.3;
else
    if ((r_block < 2.) || (r_block > 5.))
       disp('read_mars.m: r_block must be between 2.0 and 5.0!');
       return;
    end
end

fn = [ppath,filename]; disp(['filename = ',fn]);
a  = fitsread(fn);
ah = fitsinfo(fn);

xsun = cell2mat(ah.PrimaryData.Keywords(24,2));% x center
ysun = cell2mat(ah.PrimaryData.Keywords(25,2));% y center
rspix = cell2mat(ah.PrimaryData.Keywords(27,2)); %solar radius in pixels

rb = r_block*rspix;
x = [1:512] - xsun; y = [1:512] - ysun;
for i = 1:512
 for j = 1:512
    rr = sqrt(x(i)^2 + y(j)^2);
    if (rr < rb )
        a(j,i) = blocked_pix_value;
    end
 end
end

return;

        