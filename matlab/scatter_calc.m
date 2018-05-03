
%This calculates the weights for the I and Q line integrals
%  for Thomson scattering.  This was originally written for
%  ApJ 628, 1061 (2005) to produce Fig. 1.   The original file is
%  /Users/frazin/docs/Tomo/Papers/totalB_ApJ05/totalB.m


a = [1.5, 2.0,4.5] ; %impact parameters
qLimb = 0.65; %limb darkening parameter
qLimb = 0.1;
rmax = 10;

los = linspace(0,rmax,600);
wtI = zeros(length(a),length(los)); %old
wtQ = zeros(length(a),length(los)); %old

wtIn = wtI; wtQn = wtQ; %new

Ipt = wtQ;   Qpt = wtQ; % point source calculation

%test arrays
At = wtQ; Bt = At; Ct = At; Dt = At;

qi = qLimb;
for k = 1:length(a)
  r = sqrt(los.*los + a(k)^2);
  gam = asin(ones(size(los))./r);
  for i = 1:length(los)
      
      
     ri = r(i);
     gi = gam(i);
     
     %old
     ugly1 = 2*((1-qi)/(1-qi/3))*(1 - cos(gi)) + ...
           (qi/(1-qi/3))*(1 - ((cos(gi)^2)/sin(gi))*log((1+sin(gi))/cos(gi)));
     ugly2 = (2/3)*((1-qi)/(1-qi/3))*(1-cos(gi)^3) + ...
          0.25*(qi/(1-qi/3))*(1 + sin(gi)^2 - ((cos(gi)^4)/sin(gi))*log((1+sin(gi))/cos(gi)));  
     vdh = inv([2,1;2,-1])*[ugly1;ugly2];%vdh(1) = van de Hulst A parameter, (2) = B
     wtI(k,i) = ( (2-(a(k)/ri)^2)*vdh(1) + (a(k)/ri)^2*vdh(2) );
     wtQ(k,i) = ( (a(k)/ri)^2*(vdh(1) - vdh(2)) );
     
     %new  (see buildrow.c)
     rtmp  = a(k)*a(k)/(ri*ri);
     sgam = sin(gi);  cgam = cos(gi);
     vdhA = cgam * sgam * sgam;
     vdhB = (3. * sgam * sgam - 1.0 + ...
            (cgam*cgam/sgam)*(4. - 3.*cgam*cgam)*log((1.+sgam)/cgam))/8.;
     vdhC = 4./3. - cgam*(1 + cgam*cgam/3.);
     vdhD = (5. + sgam*sgam -  ... 
            (cgam*cgam/sgam)*(5. - sgam*sgam)*log((1.+sgam)/cgam))/8.;
     wtQn(k,i) = rtmp*((1.-qi)*vdhA + qi*vdhB);
     wtIn(k,i) = (1. - qi)*(2.*vdhC - vdhA*rtmp) + qi*(2.*vdhD - vdhB*rtmp);
     
     %my point source calculation
     Qpt(k,i) = .5*rtmp*rtmp;
     Ipt(k,i) = rtmp - rtmp*rtmp/2;
     
     
     At(k,i) = vdhA; Bt(k,i) = vdhB; Ct(k,i) = vdhC; Dt(k,i) = vdhD; 
     
     
  end

end % k loop over impact parameters


if 0
plot(los,wtI(1,:)/max(wtI(1,:)),'ro-',los,wtQ(1,:)/max(wtQ(1,:)),'ro:',los,wtI(2,:)/max(wtI(2,:)),'ks-',los,wtQ(2,:)/max(wtQ(2,:)),'ks:','LineWidth',3,'Markersize',10); 
   axis([0, rmax, 0, 1.05]); 
   hold on;
   set(gca,'FontSize',20);
   set(gca,'xtick',[0,2,4,6,8,10]);
   xlabel('LOS dist. (R_o)');
   ylabel('Normalized weighting');
   legend(['unpolarized, b = ',num2str(a(1)),' R_s'],['polarized, b = ',num2str(a(1)),' R_s'],['unpolarized, b = ',num2str(a(2)),' R_s'],['polarized, b = ',num2str(a(2)),' R_s']);
end
   
   
   %if print_figs
   %  eval(['print -depsc2 ',fig_dir,'tB_LOSweight']);
   %end
   
   
   