% function to compute the Green's tensor. Relevant references are:
% Aki and Richards (2002) and Haney and Nakahara (2014)

function [DC_L,DC_R] = get_sw_dc_farfield(FT,bigC,Src,Rec,max_mode,typ,pol)

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
% U_R, U_L displacements in the N,E Z (positive DOWNWARD ) coordinate
% system.
% the Green tensor is in the North-East-Down coordinate system
% To make sure the output is in meters / Newton, we have to apply a
% conversion factor of fac = 10^-12. Additionally, when we take the partial derivative of the Green's function
% we multiply by -ikr*10^-3. Total factor is then fac=10^(-15)
%% if pol = 'zrt', keep in ZRT compoennts. if pol = 'nez', convert to East North
if strcmp(pol,'rtz')==1
DC_L(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(1,length(bigC),max_mode)),...
    'U',squeeze(zeros(1,FT.nwin,max_mode)));
DC_R(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(2,length(bigC),max_mode)),...
    'U',squeeze(zeros(2,FT.nwin,max_mode)));
else
DC_L(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(2,length(bigC),max_mode)),...
    'U',squeeze(zeros(2,FT.nwin,max_mode)));
DC_R(1:Src(1).ns,1:Rec(1).nr)=struct('Uhat',squeeze(zeros(3,length(bigC),max_mode)),...
    'U',squeeze(zeros(3,FT.nwin,max_mode)));
end

fac=1E-15;

for ifreq=1:length(bigC)
    
    CC=bigC(ifreq).a;
    if isempty(CC);continue;end
    [~,ib]=unique(CC(1).zz);
    % allocate variables:
    uyi=zeros(Src(1).ns,1);uxi=uyi;uzi=uyi;duxi=uxi;duyi=uxi;duzi=uxi;
    H0L=zeros(Rec(1).nr,1);phi=H0L;C=phi;S=phi;C2=phi;S2=phi;
    ux0=H0L;uy0=ux0;uz0=ux0;
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
           
            S(ir)=sin(phi(ir));
            C(ir)=cos(phi(ir));
            % not here : Compared to Aki and Richard, the matlab convention
            % imposes that F(w) = exp(i(wt-kr)). Thus, spatial derivatives
            % of the Green's function at the source should be -ikrG.
           GL=sqrt(2/pi/CC(1).kl/Rec(ir).rg)/(8*CC(1).cl*CC(1).Ul*CC(1).Il(1))*exp(-(1i*CC(1).kl*Rec(ir).rg+pi/4));
           GR=sqrt(2/pi/CC(1).kr/Rec(ir).rg)/(8*CC(1).cr*CC(1).Ur*CC(1).Ir(1))*exp(-(1i*CC(1).kr*Rec(ir).rg+pi/4));
        % Love waves:
        toto=1i*uyi(is)*(-CC(1).kl)*( ...
            (M(1,1)-M(2,2))*S(ir)*C(ir)+M(2,1)*(S(ir)^2-C(ir)^2)) + ...
            - duyi(is)* ( M(1,3)*S(ir)-M(2,3)*C(ir)); 
        DC_L(is,ir).Uhat(ifreq) = toto*uy0(ir)*GL;

        
        toto= uxi(is)*(-CC(1).kr)*(  M(1,1)*C(ir)^2+(M(1,2)+M(2,1))*S(ir)*C(ir)...
           +M(2,2)*S(ir)^2) ...
            + 1i*duxi(is) * (  M(1,3)*C(ir) + M(2,3)*S(ir) )...
            - 1i*uzi(is) * CC(1).kr*(M(3,1)*C(ir)+M(3,2)*S(ir))+ ...
            duzi(is)*M(3,3);
        
        DC_R(is,ir).Uhat(1,ifreq) = toto*ux0(ir)*GR*exp(-1i*pi/2);
        DC_R(is,ir).Uhat(2,ifreq) = toto*uz0(ir)*GR;
        if strcmp(pol,'rtz')~=1
            toto = DC_L(is,ir).Uhat(ifreq);
            DC_L(is,ir).Uhat(1,ifreq) = -S(ir) * toto;
            DC_L(is,ir).Uhat(2,ifreq) =  C(ir) * toto ;
            
            toto = DC_R(is,ir).Uhat(1,ifreq);
            toto1 = DC_R(is,ir).Uhat(2,ifreq);
            DC_R(is,ir).Uhat(1,ifreq) = C(ir) * toto;
            DC_R(is,ir).Uhat(2,ifreq) = S(ir) * toto ;
            DC_R(is,ir).Uhat(3,ifreq) = toto1 ;
                  
        end
      
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
        S(ir)=sin(phi(ir));
        C(ir)=cos(phi(ir));
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
            GL=sqrt(2/pi/CC(1).kl(imode)/Rec(ir).rg)/(8*CC(1).cl(imode)*CC(1).Ul(imode)*CC(1).Il(1,imode))*exp(-(1i*CC(1).kl(imode)*Rec(ir).rg+pi/4));
            GR=sqrt(2/pi/CC(1).kr(imode)/Rec(ir).rg)/(8*CC(1).cr(imode)*CC(1).Ur(imode)*CC(1).Ir(1,imode))*exp(-(1i*CC(1).kr(imode)*Rec(ir).rg+pi/4));
            % Love waves:
            toto=1i*uyi(is)*(-CC(1).kl(imode))*( ...
                (M(1,1)-M(2,2))*S(ir)*C(ir)-M(2,1)*(S(ir)^2+C(ir)^2)) + ...
                - duyi(is)* ( M(1,3)*S(ir)-M(2,3)*C(ir)); 
            DC_L(is,ir).Uhat(ifreq,imode) =  toto*uy0(ir)*GL;

            toto= uxi(is)*(-CC(1).kr(imode))*(  (M(1,1)*C(ir)^2+(M(1,2)+M(2,1))*S(ir)*C(ir)...
               +M(2,2)*S(ir)^2) ...
                + 1i*duxi(is) * (  M(1,3)*C(ir)) - M(2,3)*S(ir) )...
                - 1i*uzi(is) * CC(1).kr*(M(3,1)*C(ir)+M(3,2)*S(ir))+ ...
                duzi(is)*M(3,3);
            DC_R(is,ir).Uhat(1,ifreq,imode) = -1i*toto*ux0(ir)*GR;
            DC_R(is,ir).Uhat(2,ifreq,imode) = toto*uz0(ir)*GR;
     
            
           end
       end
         if strcmp(pol,'rtz')~=1
             for is=1:SRC.ns
                 for ir=1:Rec.nr
                 phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian     
                C(ir)=cos(phi(ir));S(ir)=sin(phi(ir));
                toto = DC_L(is,ir).Uhat(:,:);
                DC_L(is,ir).Uhat(1,:,:) = -S(ir) * toto;
                DC_L(is,ir).Uhat(2,:,:) =  C(ir) * toto ;

                toto = DC_R(is,ir).Uhat(1,:,:);
                toto1 = DC_R(is,ir).Uhat(2,:,:);
                DC_R(is,ir).Uhat(1,:,:) = C(ir) * toto;
                DC_R(is,ir).Uhat(2,:,:) = S(ir) * toto ;
                DC_R(is,ir).Uhat(3,:,:) = toto1 ;
                 end
             end
                  
         end
      
    end 
    end
    if mod(ifreq/length(bigC)*100,1)<0.1
        disp([num2str(floor(ifreq/length(bigC)*100)) ' %'])
    end
end


% if requested, inverse fourier transform to provide
if strcmp(typ,'t')==1
    for is=1:Src(1).ns
        for ir=1:Rec(1).nr
        if max_mode==1
         if strcmp(pol,'rtz')==1
           % Love waves
           uu(1:length(bigC))=DC_L(is,ir).Uhat(1:length(bigC)).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_L(is,ir).U(1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
        % Rayleigh waves
          for i=1:2
           uu(1:length(bigC))=DC_R(is,ir).Uhat(i,1:length(bigC)).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_R(is,ir).U(i,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
          end

                
         else
           % Love waves
          for i=1:2
           uu(1:length(bigC))=DC_L(is,ir).Uhat(i,1:length(bigC)).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_L(is,ir).U(i,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
          end

        % Rayleigh waves
          for i=1:3
           uu(1:length(bigC))=DC_R(is,ir).Uhat(i,1:length(bigC)).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_R(is,ir).U(i,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
          end
  
        end
      else
         for imode=1:max_mode
         % Love waves
          for i=1:2
           uu(1:length(bigC))=DC_L(is,ir).Uhat(i,1:length(bigC),imode).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_L(is,ir).U(i,1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt*fac;
          end
         % Rayleigh waves
          for i=1:3
           uu(1:length(bigC))=DC_R(is,ir).Uhat(i,ifreq,imode).*Src(is).Shat(1:length(bigC));
           uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
           DC_R(is,ir).U(i,1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt*fac;
          end
         end
        end
        end
    end
end

return