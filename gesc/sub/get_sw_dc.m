% function to compute the Green's tensor. Relevant references are:
% Aki and Richards (2002) and Haney and Nakahara (2014)

function [DC_L,DC_R] = get_sw_dc(FT,bigC,Src,Rec,max_mode,typ)

%% INPUTS:
% FT : gesc format time/frequency structure
% bigC: gesc format eigen function structure
% Src: gesc source structure
%    Src(1).ns = number of sources
%    Src(:).t = time vector in seconds to go with source time function
%    Src(:).stf = source time function in time
%    Src(:).M = [M_nn, M_ne, M_nz , ... , M_zz] moment tensor
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
% U_R, U_L displacements in the N,E Z (positive downward ) coordinate
% system.
% the Green tensor is in the North-East-Down coordinate system
% To make sure the output is in meters / Newton, we have to apply a
% conversion factor of fac = 10^-12

DC_L(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(2,length(bigC),max_mode)),...
    'U',squeeze(zeros(2,FT.nwin,max_mode)));
DC_R(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(3,length(bigC),max_mode)),...
    'U',squeeze(zeros(3,FT.nwin,max_mode)));

fac=1E-12;

for ifreq=1:length(bigC)
    
    CC=bigC(ifreq).a;
    if isempty(CC);continue;end
    [~,ib]=unique(CC(1).zz);
    % allocate variables:
    uyi=zeros(Src(1).ns,1);uxi=uyi;uzi=uyi;duxi=uxi;duyi=uxi;duzi=uxi;
    H0L=zeros(Rec(1).nr,1);H1L=H0L;H2L=H0L;H3L=H0L;
    H0R=H0L;H1R=H0L;H2R=H0L;H3R=H0L;phi=H0L;C=phi;S=phi;C2=phi;S2=phi;
    ux0=H0R;uy0=ux0;uz0=ux0;
    if max_mode==1; % just fundamental
      % loop over source
       for is=1:Src(1).ns
        uyi(is) = interp1(CC(1).zz(ib),CC(1).uy(ib),Src(is).H,'spline');
        uxi(is) = interp1(CC(1).zz(ib),CC(1).ux(ib),Src(is).H,'spline');
        uzi(is) = interp1(CC(1).zz(ib),CC(1).uz(ib),Src(is).H,'spline');
        duyi(is) = -interp1(CC(1).zz(ib),CC(1).Duy(ib),Src(is).H,'spline');
        duxi(is) = -interp1(CC(1).zz(ib),CC(1).Dux(ib),Src(is).H,'spline');
        duzi(is) = -interp1(CC(1).zz(ib),CC(1).Duz(ib),Src(is).H,'spline');
       end
       % loop over receiver
       for ir=1:Rec(1).nr
         phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        % Hankel function for exact expansion of surface waves on outgoing
        % Bessel function
        H0L(ir) = besselh(0,1,CC(1).omega/CC(1).cl*Rec(ir).rg);
        H1L(ir) = besselh(1,1,CC(1).omega/CC(1).cl*Rec(ir).rg);
        H2L(ir) = besselh(2,1,CC(1).omega/CC(1).cl*Rec(ir).rg);
        H3L(ir) = besselh(3,1,CC(1).omega/CC(1).cl*Rec(ir).rg);
        H0R(ir) = besselh(0,1,CC(1).omega/CC(1).cr*Rec(ir).rg);
        H1R(ir) = besselh(1,1,CC(1).omega/CC(1).cr*Rec(ir).rg);
        H2R(ir) = besselh(2,1,CC(1).omega/CC(1).cr*Rec(ir).rg);
        H3R(ir) = besselh(3,1,CC(1).omega/CC(1).cr*Rec(ir).rg);
            % sample eigenfunctions at source and receiver depth
        if Rec(ir).H==0
         uy0(ir)=CC(1).uy(end);ux0(ir)=CC(1).ux(end);
         uz0(ir)=CC(1).uz(end);
        else
        uy0(ir) = interp1(CC(1).zz(ib),CC(1).uy(ib),Rec(ir).H,'spline');
        ux0(ir) = interp1(CC(1).zz(ib),CC(1).ux(ib),Rec(ir).H,'spline');
        uz0(ir) = interp1(CC(1).zz(ib),CC(1).uz(ib),Rec(ir).H,'spline');
        end
        S(ir)=sin(phi(ir));S2(ir)=sin(2*phi(ir));
        C(ir)=cos(phi(ir));C2(ir)=cos(2*phi(ir));
       end
       
     for is=1:Src(1).ns
         M=Src(is).M; % moment tensor in Nm
       for ir=1:Rec(1).nr
        % Love waves:
        % North component
        toto=uyi(is)*CC(1).kl/2*( ...
            (M(1,1)*C(ir)+M(1,2)*S(ir))*(H1L(ir)-C2(ir)/2*(H1L(ir)-H3L(ir))) ...
            -(M(2,1)*C(ir)+M(2,2)*S(ir))*(H1L(ir)-H3L(ir))*S2(ir)/2) ...
            + duyi(is)* ( M(1,3)/2*(H0L(ir)+H2L(ir)*C2(ir)) +M(2,3)*H2L(ir)*C(ir)*S(ir));               
        DC_L(is,ir).Uhat(1,ifreq) = 1i*toto*uy0(ir)/(8*CC(1).cl*CC(1).Ul*CC(1).Il(1));

        % East component    
        toto=uyi(is)*CC(1).kl/2*( ...
            (M(2,1)*C(ir)+M(2,2)*S(ir))*(H1L(ir)+C2(ir)/2*(H1L(ir)-H3L(ir))) ...
            -(M(1,1)*C(ir)+M(1,2)*S(ir))*(H1L(ir)-H3L(ir))*S2(ir)/2) ...
            + duyi(is)* ( M(2,3)/2*(H0L(ir)-H2L(ir)*C2(ir)) +M(1,3)*H2L(ir)*C(ir)*S(ir));      
        DC_L(is,ir).Uhat(2,ifreq) = 1i*toto*uy0(ir)/(8*CC(1).cl*CC(1).Ul*CC(1).Il(1));

        % Rayleigh waves:
        toto= uxi(is)*CC(1).kr/2*(  (M(1,1)*C(ir)+M(1,2)*S(ir))*...
            (H1R(ir)+C2(ir)/2*(H1R(ir)-H3R(ir))) +(M(2,1)*C(ir)+M(2,2)*S(ir))*...
            (H1R(ir)-H3R(ir))*S2(ir)/2) ...
            + duxi(is) * (  M(1,3)/2*(H0R(ir)-H2R(ir)*C2(ir)) - M(2,3)*H2R(ir)*S(ir)*C(ir)  )...
            - uzi(is) * CC(1).kr* (H0R(ir)-H2R(ir)) * (M(3,1)*C(ir)+M(3,2)*S(ir))*C(ir)/2 + ...
            duzi(is)*H1R(ir)*M(3,3)*C(ir);
        DC_R(is,ir).Uhat(1,ifreq) = 1i*toto*ux0(ir)/(8*CC(1).cr*CC(1).Ur*CC(1).Ir(1));

        toto= uxi(is)*CC(1).kr/2*(  (M(2,1)*C(ir)+M(2,2)*S(ir))*...
            (H1R(ir)-C2(ir)/2*(H1R(ir)-H3R(ir))) +(M(1,1)*C(ir)+M(1,2)*S(ir))*...
            (H1R(ir)-H3R(ir))*S2(ir)/2) ...
            + duxi(is) * (  M(2,3)/2*(H0R(ir)+H2R(ir)*C2(ir)) - M(1,3)*H2R(ir)*S(ir)*C(ir)  )...
            - uzi(is)*CC(1).kr * (H0R(ir)-H2R(ir)) * (M(3,1)*C(ir)+M(3,2)*S(ir))*S(ir)/2 + ...
            duzi(is)*H1R(ir)*M(3,3)*S(ir);
        DC_R(is,ir).Uhat(2,ifreq) = 1i*toto*ux0(ir)/(8*CC(1).cr*CC(1).Ur*CC(1).Ir(1));


        toto= uxi(is)*CC(1).kr/2*(H0R(ir) - H2R(ir)) * ...
            (M(1,1)*C(ir)^2 + (M(1,2)+M(2,1))*S(ir)*C(ir) +M(2,2)*S(ir)^2) ...
            - duxi(is)*H1R(ir)*(M(1,3)*C(ir)+M(2,3)*S(ir)) ...
            +uzi(is)*CC(1).kr*H1R(ir)*(M(3,1)*C(ir)+M(3,2)*S(ir)) ...
            + duzi(is)*H0R(ir)*M(3,3);
        DC_R(is,ir).Uhat(3,ifreq) = 1i*toto*uz0(ir)/(8*CC(1).cr*CC(1).Ur*CC(1).Ir(1));
        
%             crap(ifreq)= + duxi(is) ;
            
            
            %DC_R(is,ir).Uhat(3,ifreq) ;
       end
     end
    else
    for imode=1:max_mode
        
        % loop over source
       for is=1:Src(1).ns
        uyi(is) = interp1(CC(1).zz(ib),CC(1).uy(ib,imode),Src(is).H,'spline');
        uxi(is) = interp1(CC(1).zz(ib),CC(1).ux(ib,imode),Src(is).H,'spline');
        uzi(is) = interp1(CC(1).zz(ib),CC(1).uz(ib,imode),Src(is).H,'spline');
        duyi(is) = -interp1(CC(1).zz(ib),CC(1).Duy(ib,imode),Src(is).H,'spline');
        duxi(is) = -interp1(CC(1).zz(ib),CC(1).Dux(ib,imode),Src(is).H,'spline');
        duzi(is) = -interp1(CC(1).zz(ib),CC(1).Duz(ib,imode),Src(is).H,'spline');
       end
       % loop over receiver
       for ir=1:Rec(1).nr
         phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        % Hankel function for exact expansion of surface waves on outgoing
        % Bessel function
        H0L(ir) = besselh(0,1,CC(1).omega/CC(1).cl(imode)*Rec(ir).rg);
        H1L(ir) = besselh(1,1,CC(1).omega/CC(1).cl(imode)*Rec(ir).rg);
        H2L(ir) = besselh(2,1,CC(1).omega/CC(1).cl(imode)*Rec(ir).rg);
        H3L(ir) = besselh(3,1,CC(1).omega/CC(1).cl(imode)*Rec(ir).rg);
        H0R(ir) = besselh(0,1,CC(1).omega/CC(1).cr(imode)*Rec(ir).rg);
        H1R(ir) = besselh(1,1,CC(1).omega/CC(1).cr(imode)*Rec(ir).rg);
        H2R(ir) = besselh(2,1,CC(1).omega/CC(1).cr(imode)*Rec(ir).rg);
        H3R(ir) = besselh(3,1,CC(1).omega/CC(1).cr(imode)*Rec(ir).rg);
            % sample eigenfunctions at source and receiver depth
        if Rec(ir).H==0
         uy0(ir)=CC(1).uy(end,imode);ux0(ir)=CC(1).ux(end,imode);
         uz0(ir)=CC(1).uz(end,imode);
        else
         uy0(ir) = interp1(CC(1).zz(ib),CC(1).uy(ib,imode),Rec(ir).H,'spline');
         ux0(ir) = interp1(CC(1).zz(ib),CC(1).ux(ib,imode),Rec(ir).H,'spline');
         uz0(ir) = interp1(CC(1).zz(ib),CC(1).uz(ib,imode),Rec(ir).H,'spline');
        end
       end
       % Build displacement spectra
       for is=1:Src(1).ns
         M=Src(is).M; % moment tensor in Nm
           for ir=1:Rec(1).nr
            % Love waves:
            % North component
            toto=uyi(is)*CC(1).kl(imode)/2*( ...
                (M(1,1)*C(ir)+M(1,2)*S(ir))*(H1L(ir)-C2(ir)/2*(H1L(ir)-H3L(ir))) ...
                -(M(2,1)*C(ir)+M(2,2)*S(ir))*(H1L(ir)-H3L(ir))*S2(ir)/2) ...
                + duyi(is)* ( M(1,3)/2*(H0L(ir)+H2L(ir)*C2(ir)) +M(2,3)*H2L(ir)*C(ir)*S(ir));               
            DC_L(is,ir).Uhat(1,ifreq) = 1i*toto*uy0(ir)/(8*CC(1).cl(imode)*CC(1).Ul(imode)*CC(1).Il(1,imode));

            % East component    
            toto=uyi(is)*CC(1).kl(imode)/2*( ...
                (M(2,1)*C(ir)+M(2,2)*S(ir))*(H1L(ir)+C2(ir)/2*(H1L(ir)-H3L(ir))) ...
                -(M(1,1)*C(ir)+M(1,2)*S(ir))*(H1L(ir)-H3L(ir))*S2(ir)/2) ...
                + duyi(is)* ( M(2,3)/2*(H0L(ir)-H2L(ir)*C2(ir)) +M(1,3)*H2L(ir)*C(ir)*S(ir));      
            DC_L(is,ir).Uhat(2,ifreq) = 1i*toto*uy0(ir)/(8*CC(1).cl(imode)*CC(1).Ul(imode)*CC(1).Il(1,imode));
            
            % Rayleigh waves:
        toto= uxi(is)*CC(1).kr(imode)/2*(  (M(1,1)*C(ir)+M(1,2)*S(ir))*...
            (H1R(ir)+C2(ir)/2*(H1R(ir)-H3R(ir))) +(M(2,1)*C(ir)+M(2,2)*S(ir))*...
            (H1R(ir)-H3R(ir))*S2(ir)/2) ...
            + duxi(is) * (  M(1,3)/2*(H0R(ir)-H2R(ir)*C2(ir)) - M(2,3)*H2R(ir)*S(ir)*C(ir)  )...
            - uzi(is) * CC(1).kr(imode)* (H0R(ir)-H2R(ir)) * (M(3,1)*C(ir)+M(3,2)*S(ir))*C(ir)/2 + ...
            duzi(is)*H1R(ir)*M(3,3)*C(ir);
            DC_R(is,ir).Uhat(1,ifreq) = 1i*toto*ux0(ir)/(8*CC(1).cr(imode)*CC(1).Ur(imode)*CC(1).Ir(1,imode));
            
        toto= uxi(is)*CC(1).kr(imode)/2*(  (M(2,1)*C(ir)+M(2,2)*S(ir))*...
            (H1R(ir)-C2(ir)/2*(H1R(ir)-H3R(ir))) +(M(1,1)*C(ir)+M(1,2)*S(ir))*...
            (H1R(ir)-H3R(ir))*S2(ir)/2) ...
            + duxi(is) * (  M(2,3)/2*(H0R(ir)+H2R(ir)*C2(ir)) - M(1,3)*H2R(ir)*S(ir)*C(ir)  )...
            - uzi(is)*CC(1).kr(imode) * (H0R(ir)-H2R(ir)) * (M(3,1)*C(ir)+M(3,2)*S(ir))*S(ir)/2 + ...
            duzi(is)*H1R(ir)*M(3,3)*S(ir);
            DC_R(is,ir).Uhat(2,ifreq) = 1i*toto*ux0(ir)/(8*CC(1).cr(imode)*CC(1).Ur(imode)*CC(1).Ir(1,imode));
            
            
        toto= uxi(is)*CC(1).kr(imode)/2*(H0R(ir) - H2R(ir)) * ...
            (M(1,1)*C(ir)^2 + (M(1,2)+M(2,1))*S(ir)*C(ir) +M(2,2)*S(ir)^2) ...
            - duxi(is)*H1R(ir)*(M(1,3)*C(ir)+M(2,3)*S(ir)) ...
            +uzi(is)*CC(1).kr(imode)*H1R(ir)*(M(3,1)*C(ir)+M(3,2)*S(ir)) ...
            + duzi(is)*H0R(ir)*M(3,3);
            DC_R(is,ir).Uhat(3,ifreq) = 1i*toto*uz0(ir)/(8*CC(1).cr(imode)*CC(1).Ur(imode)*CC(1).Ir(1,imode));
            
            
           end
       end
    end
    end
%     disp([num2str(ifreq/length(bigC)*100) ' %'])
end

% if requested, inverse fourier transform to provide
if strcmp(typ,'t')==1
    
tt=linspace(0,FT.twin,FT.nwin);

    for is=1:Src(1).ns
%         ik=find(Src(is).t(end)>=tt);
%     SS = interp1(Src(is).t,Src(is).stf,tt(ik),'spline');
%     SS(ik(end)+1:FT.nwin)=0;
%     Src(is).Shat=fft(SS);
      for ir=1:Rec(1).nr
        if max_mode==1
         % Love waves
          for i=1:2
           uu(1:length(bigC))=conj(DC_L(is,ir).Uhat(i,1:length(bigC))).*(Src(is).Shat(1:length(bigC)));
           uu(length(bigC)+1:FT.nwin)=0;
           DC_L(is,ir).U(i,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt;
          end

        % Rayleigh waves
          for i=1:3
           uu(1:length(bigC))=conj(DC_R(is,ir).Uhat(i,1:length(bigC))).*(Src(is).Shat(1:length(bigC)));
           uu(length(bigC)+1:FT.nwin)=0;
           DC_R(is,ir).U(i,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt;
          end

        else
         for imode=1:max_mode
         % Love waves
          for i=1:2
           uu(1:length(bigC))=conj(DC_L(is,ir).Uhat(i,1:length(bigC),imode)).*(Src(is).Shat(1:length(bigC)));
           uu(length(bigC)+1:FT.nwin)=0;
           DC_L(is,ir).U(i,1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt;
          end
         % Rayleigh waves
          for i=1:3
           uu(1:length(bigC))=conj(DC_R(is,ir).Uhat(i,ifreq,imode)).*(Src(is).Shat(1:length(bigC)));
           uu(length(bigC)+1:FT.nwin)=0;
           DC_R(is,ir).U(i,1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt;
          end
         end
       end
        DC_R(is,ir).U=DC_R(is,ir).U*fac;
        DC_L(is,ir).U=DC_L(is,ir).U*fac;
      end
    end
end

% 
% 
% figure(12)
% loglog(FT.omega(1:length(bigC)),abs(crap));
return