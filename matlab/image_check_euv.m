function  medianlist = image_check_euv(listfile,ppath)
%function medianlist = image_check_euv(listfile,ppath)
%
% IF THERE IS A PROBLEM READING FILES SEE THE LINE:
%     crap = crap(2:length(crap)); and adjust as needed
%
%this reads a list file for euv tomography and displays the images
%   
%medianlist - median of each image
%listfile - file containing list of filenames.  first line is the number of
%    files
%ppath - directory where filename resides.  include final slash



if (nargin ~= 2)
    disp('image_check_euv.m: Wrong number of input arguments!');
    return;
end

fn = [ppath,listfile]; %disp(['filename = ',fn]);

fid = fopen(fn,'r');
if (fid < 0)
    error(['cannot open file: ',fn]);
end

t = textscan(fid,'%4s',1);

nfiles = str2num(cell2mat(t{1}(1))); % the first cell is the number of files
disp(['image_check_euv.m: there are ',num2str(nfiles),' files.']);

medianlist = zeros(1,nfiles);% keep track of image median values
fnn = cell(1,nfiles);% keep track of filenames
carlong = cell(1,nfiles); %carrington longitudes from fits files
carlat  = carlong; % latitudes

figure;
for k = 1:nfiles
    t = textscan(fid,'%80s',1);
    crap = cell2mat(t{1}(1));
    if (k == -1)
        crap = crap(2:length(crap));
    end
    fnn{k} = crap ;
    [a,h] = read_euv_image(fnn{k},ppath);
    
    hdr_ind =  find(strcmp(h.PrimaryData.Keywords(:,1),'CRLN_OBS')==1);
    carlong(k) = h.PrimaryData.Keywords(hdr_ind,2);
    hdr_ind =  find(strcmp(h.PrimaryData.Keywords(:,1),'CRLT_OBS')==1);
    carlat(k)  = h.PrimaryData.Keywords(hdr_ind,2);
    
    medianlist(k) = median(a(:));
    a( find(a < 0) ) = 0;
    imagesc(sqrt(a/medianlist(k)));
    set(gca,'YDir','normal'); colorbar;
    title(['long = ',carlong(k),' median = ',num2str(medianlist(k))]);
    xlabel(['file ',num2str(k),': ',fnn{k}]);
    pause;
end
fclose(fid);
clear crap
