function [chisq,y_orig,y_synth] = chisq_stat(lambda,fn_orig,fn_synth,dir)
%function [chisq,y_orig,y_synth] = chisq_stat(fn_orig,fn_synth,dir)
%This calculates a chi-squared statistic for a tomographic y vectors.
%The inputs are the filenames of the y vector containing the original
% data and the synthetic data.  So far, only poisson stats are included
%
% chisqu is the chi-squared value
% y_orig, y_synth are the original and synthetic y vectors
% lambda = '171', '195' or '284'
% fn_orig is the filename of the original y vector
% fn_synth                       synthetic
% dir is the directory containing the vectors.  include final /

ndf = 13*60*120;% number of degrees of freedom
tbin = 3; % time binning factor
pbin = 4; % 1D pixel binning factor 

fno = [dir, fn_orig];
fns = [dir, fn_synth]; 

fid = fopen(fno,'rb');
if (fid == -1)
    disp(['cannot open file: ', fno]);
    return;
else
    y_orig = fread(fid,inf,'float32');
    fclose(fid);
end


fid = fopen(fns,'rb');
if (fid == -1)
    disp(['cannot open file: ', fns]);
    return;
else
    y_synth = fread(fid,inf,'float32');
    fclose(fid);
end

if (length(y_synth) ~= length(y_orig))
   disp('file lengths do not match!');
   disp(['length(yo) = ',num2str(length(yo)),'length(ys) = ',num2str(length(ys))]);
   return;
end

if (lambda == '171')
    const = 0.675;
elseif (lambda == '195')
    const = 0.148;
elseif (lambda == '284')
    const = 0.0754;
else
    disp(['lambda not set properly']);
    return;
end

sigmasq = (y_orig*const/(pbin*pbin*tbin));

chisq = sum(((y_orig - y_synth).^2)./sigmasq);
chisq = chisq/(1 + length(y_orig) - ndf);





