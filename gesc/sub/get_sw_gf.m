% function to compute the Green's tensor. Revelant references are:
% Aki and Richards (2002) and Haney and Nakahara (2014)

function [G_L,G_R] = get_sw_gf(FT,bigC,Src,Rec,max_mode,typ)

%% INPUTS:
% FT : gesc format time/frequency structure
% bigC: gesc format eigen function structure
% Src: gesc source structure
%    Src(1).ns = number of sources
%    Src(:).t = time vector in seconds to go with source time function
%    Src(:).stf = source time function in time
%    Src(:).F = [F_n F_e F_z] force vector in Newton
%    Src(:).HS: source depth * in km*
% Rec: gesc structure for receiver:
%    Rec(1).nr = number of receiver
%    Rec(:).rg: distance from source *im km*
%    Rec(:).az: azimuth from source *im degrees*
%    Rec(:).H: depth of receiver *im km*
% max_mode: number of modes
% typ: 'f', output in Fourier Domain
%      't', output in time domain
%% OUTPUTS:
% G_L and G_R are structure of size (Ns,Nr) and that contains the Green's
% tensor components (2x2 for Love waves, and 3x3 for Rayleigh waves)
% in (m/N)
% the Green tensor is in the North-East-Down coordinate system
% To make sure the output is in meters / Newton, we have to apply a
% conversion factor of fac = 10^-12


fac=1E-12;

for ifreq=1:length(bigC)
    
    C=bigC(ifreq).a;
    if isempty(C);continue;end
    [~,ib]=unique(C(1).zz);
    if max_mode==1; % just fundamental
      % loop over source
       for is=1:Src(1).ns
        uyi(is) = interp1(C(1).zz(ib),C(1).uy(ib),Src(is).H,'spline');
        uxi(is) = interp1(C(1).zz(ib),C(1).ux(ib),Src(is).H,'spline');
        uzi(is) = interp1(C(1).zz(ib),C(1).uz(ib),Src(is).H,'spline');
       end
       % loop over receiver
       for ir=1:Rec(1).nr
         phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        % Hankel function for exact expansion of surface waves on outgoing
        % Bessel function
        H01L(ir) = besselh(0,1,C(1).omega/C(1).cl*Rec(ir).rg);
        H21L(ir) = besselh(2,1,C(1).omega/C(1).cl*Rec(ir).rg);
        H01R(ir) = besselh(0,1,C(1).omega/C(1).cr*Rec(ir).rg);
        H11R(ir) = besselh(1,1,C(1).omega/C(1).cr*Rec(ir).rg);
        H21R(ir) = besselh(2,1,C(1).omega/C(1).cr*Rec(ir).rg);
            % sample eigenfunctions at source and receiver depth
        if Rec(ir).H==0
         uy0(ir)=C(1).uy(end);ux0(ir)=C(1).ux(end);
         uz0(ir)=C(1).uz(end);
        else
        uy0(ir) = interp1(C(1).zz(ib),C(1).uy(ib),Rec(ir).H,'spline');
        ux0(ir) = interp1(C(1).zz(ib),C(1).ux(ib),Rec(ir).H,'spline');
        uz0(ir) = interp1(C(1).zz(ib),C(1).uz(ib),Rec(ir).H,'spline');
        end
       end
       % build Green's function
       for is=1:Src(1).ns
           for ir=1:Rec(1).nr
            % Love waves:
            G(1,1) = 1/2*(H01L(ir)+H21L(ir)*cos(2*phi(ir)));
            G(1,2) = H21L(ir)*sin(phi(ir))*cos(phi(ir));
            G(2,1) = H21L(ir)*sin(phi(ir))*cos(phi(ir));
            G(2,2) = 1/2*(H01L(ir)-H21L(ir)*cos(2*phi(ir)));
            
            G_L(is,ir).Ghat(1:2,1:2,ifreq)=1i*uy0(ir)*uyi(is)/...
                (8*C(1).cl*C(1).Ul*C(1).Il(1))*G;
            % Rayleigh waves:
            Gr(1,1) = ux0(ir)*uxi(is)/2*(H01R(ir)-H21R(ir)*cos(2*phi(ir)));
            Gr(1,2) = -ux0(ir)*uxi(is)*H21R(ir)*sin(phi(ir))*cos(phi(ir));
            Gr(1,3) = ux0(ir)*uzi(is)*H11R(ir)*cos(phi(ir));
            Gr(2,1) = -ux0(ir)*uxi(is)*H21R(ir)*sin(phi(ir))*cos(phi(ir));
            Gr(2,2) = ux0(ir)*uxi(is)/2*(H01R(ir)+H21R(ir)*cos(2*phi(ir)));
            Gr(2,3) = ux0(ir)*uzi(is)*H11R(ir)*sin(phi(ir));
            Gr(3,1) = -uz0(ir)*uxi(is)*H11R(ir)*cos(phi(ir));
            Gr(3,2) = -uz0(ir)*uxi(is)*H11R(ir)*sin(phi(ir));
            Gr(3,3) = uz0(ir)*uzi(is)*H01R(ir);
            G_R(is,ir).Ghat(1:3,1:3,ifreq) = Gr*1i/(8*C(1).cr*C(1).Ur*C(1).Ir(1));
           end
       end
     
    else
        for imode=1:max_mode
        
        % loop over source
       for is=1:Src(1).ns
        uyi(is) = interp1(C(1).zz(ib),C(1).uy(ib,imode),Src(is).H,'spline');
        uxi(is) = interp1(C(1).zz(ib),C(1).ux(ib,imode),Src(is).H,'spline');
        uzi(is) = interp1(C(1).zz(ib),C(1).uz(ib,imode),Src(is).H,'spline');
       end
       % loop over receiver
       for ir=1:Rec(1).nr
         phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        % Hankel function for exact expansion of surface waves on outgoing
        % Bessel function
        H01L(ir) = besselh(0,1,C(1).omega/C(1).cl(imode)*Rec(ir).rg);
        H21L(ir) = besselh(2,1,C(1).omega/C(1).cl(imode)*Rec(ir).rg);
        H01R(ir) = besselh(0,1,C(1).omega/C(1).cr(imode)*Rec(ir).rg);
        H11R(ir) = besselh(1,1,C(1).omega/C(1).cr(imode)*Rec(ir).rg);
        H21R(ir) = besselh(2,1,C(1).omega/C(1).cr(imode)*Rec(ir).rg);
            % sample eigenfunctions at source and receiver depth
        if Rec(ir).H==0
         uy0(ir)=C(1).uy(end,imode);ux0(ir)=C(1).ux(end,imode);
         uz0(ir)=C(1).uz(end,imode);
        else
        uy0(ir) = interp1(C(1).zz(ib),C(1).uy(ib,imode),Rec(ir).H,'spline');
        ux0(ir) = interp1(C(1).zz(ib),C(1).ux(ib,imode),Rec(ir).H,'spline');
        uz0(ir) = interp1(C(1).zz(ib),C(1).uz(ib,imode),Rec(ir).H,'spline');
        end
       end
       % build Green's function
       for is=1:Src(1).ns
           for ir=1:Rec(1).nr
            % Love waves:
            G(1,1) = 1/2*(H01L(ir)+H21L(ir)*cos(2*phi(ir)));
            G(1,2) = H21L(ir)*sin(phi(ir))*cos(phi(ir));
            G(2,1) = H21L(ir)*sin(phi(ir))*cos(phi(ir));
            G(2,2) = 1/2*(H01L(ir)-H21L(ir)*cos(2*phi(ir)));
            G_L(is,ir).Ghat(1:2,1:2,ifreq,imode)=1i*uy0(ir)*uyi(is)/...
                (8*C(1).cl(imode)*C(1).Ul(imode)*C(1).Il(1,imode))*G;

            % Rayleigh waves:
            Gr(1,1) = ux0(ir)*uyi(is)/2*(H01R(ir)-H21R(ir)*cos(2*phi(ir)));
            Gr(1,2) = -ux0(ir)*uxi(is)*H21R(ir)*sin(phi(ir))*cos(phi(ir));
            Gr(1,3) = ux0(ir)*uzi(is)*H11R(ir)*cos(phi(ir));
            Gr(2,1) = -ux0(ir)*uxi(is)*H21R(ir)*sin(phi(ir))*cos(phi(ir));
            Gr(2,2) = ux0(ir)*uxi(is)/2*(H01R(ir)+H21R(ir)*cos(2*phi(ir)));
            Gr(2,3) = ux0(ir)*uzi(is)*H11R(ir)*sin(phi(ir));
            Gr(3,1) = -uz0(ir)*uxi(is)*H11R(ir)*cos(phi(ir));
            Gr(3,2) = -uz0(ir)*uxi(is)*H11R(ir)*sin(phi(ir));
            Gr(3,3) = uz0(ir)*uzi(is)*H01R(ir);
            G_R(is,ir).Ghat(1:3,1:3,ifreq,imode) = Gr*1i/(8*C(1).cr(imode)*C(1).Ur(imode)*C(1).Ir(1,imode));
           end
       end
    end
    end
    disp([num2str(ifreq/length(bigC)*100) ' %'])
end


% if requested, inverse fourier transform to provide
if strcmp(typ,'t')==1
    for is=1:Src(1).ns
        for ir=1:Rec(1).nr
         if max_mode==1
                % Love waves
       for ii=1:2
           for jj=1:2
               uu(1:length(bigC))=G_L(is,ir).Ghat(ii,jj,1:length(bigC));
               uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
               G_L(is,ir).G(ii,jj,1:FT.nwin)=ifft(conj(uu),'symmetric')/FT.dt;
           end
       end
        % Rayleigh waves
       for ii=1:3
           for jj=1:3
               uu(1:length(bigC))=G_R(is,ir).Ghat(ii,jj,1:length(bigC));
               uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
               G_R(is,ir).G(ii,jj,1:FT.nwin)=ifft(conj(uu),'symmetric')/FT.dt;
           end
       end
                
       
            else
        for imode=1:max_mode
        % Love waves
       for ii=1:2
           for jj=1:2
               uu(1:length(bigC))=G_L(is,ir).Ghat(ii,jj,1:length(bigC),imode);
               uu(length(bigC)+1:FT.nwin)=0;
               G_L(is,ir).G(ii,jj,1:FT.nwin,imode)=ifft(conj(uu),'symmetric')/FT.dt;
           end
       end
        % Rayleigh waves
       for ii=1:3
           for jj=1:3
               uu(1:length(bigC))=G_R(is,ir).Ghat(ii,jj,ifreq,imode);
               uu(length(bigC)+1:FT.nwin)=0;
               G_R(is,ir).G(ii,jj,1:FT.nwin,imode)=ifft(conj(uu),'symmetric')/FT.dt;
           end
       end
        end
            end
            G_R(is,ir).G=G_R(is,ir).G*fac;
            G_L(is,ir).G=G_L(is,ir).G*fac;
        end
    end
end


return