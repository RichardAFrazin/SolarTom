function imchk = image_check_c2(listfile,ppath,r_block)
%function imchk = image_check_c2(listfile,ppath,r_block)
%
% IF THERE IS A PROBLEM READING FILES SEE THE LINE:
%     crap = crap(2:length(crap)); and adjust as needed
%
%this reads a list file for tomography and displays the images
%   This reads the Marseilles LASCO pB images.
%imchk - output
%
%ppath - directory where filename resides.  include final slash
%r_block - optional argument for the blocking radius in units of
%   Rsun where Rsun is taken to be the instaneously correct value. Default
%   value is 2.3
%see also: read_mars.m
%


if ((nargin < 2) || (nargin > 3))
    disp('image_check: Wrong number of input arguments!');
    return;
end
if nargin == 2
    r_block = 2.3;
else
    if ((r_block < 2.) || (r_block > 5.))
       disp('image_check: r_block must be between 2.0 and 5.0!');
       return;
    end
end

fn = [ppath,listfile]; %disp(['filename = ',fn]);

fid = fopen(fn,'r');
t = textscan(fid,'%4s',1);

nfiles = str2num(cell2mat(t{1}(1))); % the first cell is the number of files
disp(['image_check.m: there are ',num2str(nfiles),' files.']);

imchk = zeros(512,512,nfiles);
filter = cell(nfiles);
fnn = cell(1,nfiles);
for k = 1:nfiles
    t = textscan(fid,'%35s',1);
    crap = cell2mat(t{1}(1));
    if (k == -1)
        crap = crap(2:length(crap));
    end
    fnn{k} = crap ;
    [a,fltr,h] = read_mars_nrl(fnn{k},ppath,r_block);
    imchk(:,:,k) = a;
    filter{k} = fltr;
end
fclose(fid);
clear crap

cmin = min(min(min(a))); cmax = max(max(max(a)));

figure;
for k = 1:nfiles
   imagesc(sqrt(imchk(:,:,k))/sqrt(cmax)); caxis([sqrt(cmin), 1]); set(gca,'YDir','normal');
   title([num2str(k),': ',fnn{k},' ',filter{k}]); colorbar;
   if ~strcmpi(filter{k},'Orange')
      disp(['WARNING - bad filter: ',filter{k},' ',fnn{k}]);
   end
   pause;
end
    
end


