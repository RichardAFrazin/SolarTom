%
function r = snr4D(x,x0,iradmin,iradmax,nrad,ntheta,nphi,ntime);
    
    x=reshape(x,nrad,ntheta,nphi,ntime);
    x0=reshape(x0,nrad,ntheta,nphi,ntime);
    x0cut = x0(iradmin:iradmax,:,:,:);
    xcut = x(iradmin:iradmax,:,:,:);
    r = mean(xcut(:))/std(xcut(:)-x0cut(:));
    % r = 20 * log10(r) % en dB
end
