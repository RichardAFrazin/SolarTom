function [a,filter,ah] = read_mars_nrl(filename,ppath,r_block)

%function read_mars_nrl(filename,r_block)
%   This reads the Marseilles and NRL C2 images.
%a - 512x512 output image with blocking
%filter - the 'FILTER' string from the header structure
%ah - fits header structure
%filename  - duh
%ppath - directory where filename resides.  include final slash
%r_block - optional argument for the blocking radius in units of
%   Rsun where Rsun is taken to be the instaneously correct value. Default
%   value is 2.3

%blocked_pix_value = 0.1;

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

filteri =  find(strcmp(ah.PrimaryData.Keywords(:,1),'FILTER')==1);
filter = cell2mat(ah.PrimaryData.Keywords(filteri,2));

NRLx = find(strcmp(ah.PrimaryData.Keywords(:,1),'CRPIX1')==1);
NRLy = find(strcmp(ah.PrimaryData.Keywords(:,1),'CRPIX2')==1);

MARx = find(strcmp(ah.PrimaryData.Keywords(:,1),'XSUN')==1);
MARy = find(strcmp(ah.PrimaryData.Keywords(:,1),'YSUN')==1);

RSUN = find(strcmp(ah.PrimaryData.Keywords(:,1),'RSUN_PIX')==1);

if(isempty(MARx) || isempty(MARy)) 
    disp('Detected NRL File....');
    xsun =  cell2mat(ah.PrimaryData.Keywords(NRLx,2));% x center
    ysun =  cell2mat(ah.PrimaryData.Keywords(NRLy,2));% y center
else
    disp('Detected Marseilles File....');
    xsun =  cell2mat(ah.PrimaryData.Keywords(MARx,2));% x center
    ysun =  cell2mat(ah.PrimaryData.Keywords(MARy,2));% y center
end

if(isempty(RSUN))
    disp('Could not find solar radius... picking value');
    rspix = 40;
else
    rspix = cell2mat(ah.PrimaryData.Keywords(RSUN,2)); %solar radius in pixels
end

blocked_pix_value = min(min(a));


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

        