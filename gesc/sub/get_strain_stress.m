function [EE, SS] = get_strain_stress(MED,FT,bigC,R,T,Z,H,or);

%% the tensors are in the radial-transverse-down orientation
% EE is strain tensor of dimension: 3 x 3 x # of time samples 
% SS is stress tensor of dimension: 3 x 3 x # of time samples (in Pascal)
% MED is the medium structure provided by GESC.
% FT is the Fourier-Time structure provvided by GESC.
% Rec is the receiver structures in the main code
% R, T, Z, are the radial transverse and vertical components of ground
% motion with the dimension # of receivers x # of time samples
% H is the depth to wich evaluate the strain and stress
% or is 'rtz' for the cynlindrical coordinate of 'nez' for the
% North-East-Down coordinate.

% EE=zeros(3,3,length(Z(:,1)));SS=EE;
EE=struct('RR',zeros(length(Z(:,1)),FT.nwin),'TT',zeros(length(Z(:,1)),FT.nwin),'ZZ',zeros(length(Z(:,1)),FT.nwin),...
    'RT',zeros(length(Z(:,1)),FT.nwin),'TZ',zeros(length(Z(:,1)),FT.nwin),'RZ',zeros(length(Z(:,1)),FT.nwin));
SS=EE;
alpha = interp1(MED(1).z,MED(1).alpha,H,'linear')*1E3;
beta = interp1(MED(1).z,MED(1).beta, H,'linear')*1E3;
rho = interp1(MED(1).z,MED(1).rho,H,'linear')*1E3;
mu = rho*beta^2;
lambda  = rho.*alpha.^2 - 2*mu;

for ir=1:length(Z(:,1))
    Rh(ir,:)=fft(R(ir,:));
    Th(ir,:)=fft(T(ir,:));
    Zh(ir,:)=fft(Z(ir,:));
end
for ifreq=2:length(bigC)
    
    C=bigC(ifreq).a;
    if isempty(C);continue;end
    [~,ib]=unique(C(1).zz);
    uyi(ifreq) = interp1(C(1).zz(ib),C(1).uy(ib),H,'spline');
    uxi(ifreq) = interp1(C(1).zz(ib),C(1).ux(ib),H,'spline');
    uzi(ifreq) = interp1(C(1).zz(ib),C(1).uz(ib),H,'spline');
    duxi(ifreq) = -interp1(C(1).zz(ib),C(1).Dux(ib,1),H,'spline');
    duzi(ifreq) = -interp1(C(1).zz(ib),C(1).Duz(ib,1),H,'spline');
    duyi(ifreq) = -interp1(C(1).zz(ib),C(1).Duy(ib,1),H,'spline');
    uy0(ifreq)=C(1).uy(end);ux0(ifreq)=C(1).ux(end);
    uz0(ifreq)=C(1).uz(end);
    kr(ifreq) = C(1).kr(1);
    kl(ifreq) = C(1).kl(1);
end
II=(2:length(bigC));
for ir=1:length(Z(:,1))
    crap(II)=1i*kr(II).*uxi(II)./ux0(II).*Rh(ir,II);
    crap(II(end)+1:FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.RR(ir,1:FT.nwin)= ifft(crap,'symmetric')*1E-3;
    
%     crap(II)=1i*kr.*uxi./ux0.*Rh(ir,II);
%     crap(II(end)+1+FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.TT(ir,1:FT.nwin)= zeros(FT.nwin,1);%ifft(crap,'symmetric');
    
    crap(II)= duzi(II)./uz0(II).*Zh(ir,II);
    crap(II(end)+1:FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.ZZ(ir,1:FT.nwin)= ifft(crap,'symmetric')*1E-3;
    
    
    crap(II)= 1/2*1i*kl(II).*Th(ir,II);
    crap(II(end)+1:FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.RT(ir,1:FT.nwin)= ifft(crap,'symmetric')*1E-3;

    crap(II)=  1/2*duyi(II)./uy0(II).*Th(ir,II);
    crap(II(end)+1:FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.TZ(ir,1:FT.nwin)= ifft(crap,'symmetric')*1E-3;

    
    crap(II)= 1/2*duxi(II)./ux0(II).*R(ir,II) + 1/2*1i*kr(II).*Zh(ir,II);
    crap(II(end)+1:FT.nwin)=0;crap(1)=0;crap(II)=crap(II)-mean(crap(II));    
    EE.RZ(ir,1:FT.nwin)= ifft(crap,'symmetric')*1E-3;  
end
SS.RR = lambda*(EE.RR+EE.TT+EE.ZZ) + 2*mu*EE.RR;
SS.TT = lambda*(EE.RR+EE.TT+EE.ZZ) + 2*mu*EE.TT;
SS.ZZ = lambda*(EE.RR+EE.TT+EE.ZZ) + 2*mu*EE.ZZ;
SS.RT = 2*mu*EE.RT;
SS.RZ = 2*mu*EE.RZ;
SS.TZ = 2*mu*EE.TZ;

if strcmp(or,'nez')==1
   
        SS.NN = cosa*cosb*T(1,:)-cosa*sinb*T(2,:)+sina*cosb*T(4,:)-sina*sinb*T(5,:);
        SS.NE = cosa*sinb*T(1,:)+cosa*cosb*T(2,:)-sina*sinb*T(4,:)-sina*cosb*T(5,:);
        TOTO1{4}(npair,:) = sina*cosb*T(1,:)-sina*sinb*T(2,:)+cosa*cosb*T(4,:)-cosa*sinb*T(5,:);
        TOTO1{5}(npair,:) = sina*sinb*T(1,:)+sina*cosb*T(2,:)+cosa*sinb*T(4,:)+cosa*cosb*T(5,:);
        TOTO1{3}(npair,:) = cosa*T(3,:)-sina*T(6,:);
        TOTO1{6}(npair,:) = sina*T(3,:)+cosa*T(6,:);
        TOTO1{7}(npair,:) = cosb*T(7,:)-sinb*T(8,:);
        TOTO1{8}(npair,:) = sinb*T(7,:)+cosb*T(8,:);
        TOTO1{9}(npair,:) = T(9,:); 
    
end
end