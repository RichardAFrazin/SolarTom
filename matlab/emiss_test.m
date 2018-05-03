pth = '/Users/frazin/ML/Solar/Data/Alberto_Mar07/';
fid = fopen([pth,'responses_800.dat'],'rb');
a = fread(fid,inf,'float32','ieee-le');
fclose(fid);
T = a(1:200); r171 = a(201:400); r195 = a(401:600); r284 = a(601:800);
r171 = r171/3.1e-31; r195 = r195/3.1e-31; r284 = r284/3.1e-33;
T = T/1.e6;

%Gaussian parameters of 'true' signal
%Tc = 1.1; sc = .3; ar = 1;
%xt =   (ar/2/pi/sc)*exp(- ((T - Tc).^2)/(2.*sc^2) );
xt = zeros(length(T));
xt(95:115) = 2.5;


%emissvity
e171 = r171'*xt; e195 = r195'*xt; e284 = r284'*xt;

%set up grid
Tcgrid = logspace(-.5,.6,25); lTc = length(Tcgrid);
scgrid = logspace(-1,log10(2),26); lsc = length(scgrid);
argrid = logspace(-.5,1.5,20); lar = length(argrid);

%cost function
phi = zeros(lTc, lsc, lar);
for j = 1:lTc
  for k = 1:lsc
    for l = 1:lar
        xTc = Tcgrid(j); xsc = scgrid(k); xar = argrid(l);
        f = (xar/2/pi/xsc)*exp( - ((T - xTc).^2)/(2.*xsc^2) );
        phi(j,k,l) = (e171 - r171'*f)^2  + ...
            (e195 - r195'*f)^2 + (e284 - r284'*f)^2;
    end
  end
end

figure(1);
for l= 1:lar
    crap = phi(:,:,l);
    %fcrap = find(crap > 30);
    %crap(fcrap) = 30;
    imagesc(log10(crap)); colorbar; 
    title([num2str(argrid(l)),' ', num2str(min(min(crap)))]);
    pause;
end

