%
function r = snr(x,x0,iradmax);
    
    x=reshape(x,60,60,120);
    x0=reshape(x0,60,60,120);
    x0cut = x0(1:iradmax,:,:);
    xcut = x(1:iradmax,:,:);
    r = mean(xcut(:))/std(xcut(:)-x0cut(:));
    % r = 20 * log10(r) % en dB
end
