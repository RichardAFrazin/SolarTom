function [a, ah] = read_euv_image(filename,ppath);
%function [a, ah] = read_euv_image(filename,ppath);
%
%a  - image data
%ah - header structure
%filename - self-explanatory
%ppath - path to file, use final slash

fn = [ppath,filename]; disp(['filename = ',fn]);
a  = fitsread(fn);
ah = fitsinfo(fn);

