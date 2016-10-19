% function to compute the Green's tensor. Revelant references are:
% Aki and Richards (2002) and Haney and Nakahara (2014)

function [G_L,G_R] = get_sw_gf_farfield(FT,bigC,Src,Rec,max_mode,typ,pol)

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
% tensor components (TT for Love waves, and RR,RZ,ZR,ZZ for Rayleigh waves)
% in (m/N)
% the Green tensor is in the North-East-Down coordinate system
% To make sure the output is in meters / Newton, we have to apply a
% conversion factor of fac = 10^-12
fac=1E-12;

%% NOTES on MATLAB FFT
% matlab's convention for fft is:
%F(w) = int_{-\infty}^{\infty} f(t) exp(-iwt) dt
%f(t) = int_{-\infty}^{\infty} F(w) exp(iwt) dt
% To have a forward propagating wave , we need to use the convention
% F(w) = exp(-ikr).

for ifreq=1:length(bigC)
    C=bigC(ifreq).a;
    if isempty(C);continue;end
    [~,ib]=unique(C(1).zz);
    if max_mode==1; % just fundamental
      %% loop over source
       for is=1:Src(1).ns
           if Src(is).H==0
               uyi(is)=C(1).uy(end);
               uxi(is)=C(1).ux(end);
               uzi(is)=C(1).uz(end);
           else
               uyi(is) = interp1(C(1).zz(ib),C(1).uy(ib),Src(is).H,'linear');
               uxi(is) = interp1(C(1).zz(ib),C(1).ux(ib),Src(is).H,'linear');
               uzi(is) = interp1(C(1).zz(ib),C(1).uz(ib),Src(is).H,'linear');
           end
       end
       %% loop over receiver
       for ir=1:Rec(1).nr
        phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        J0R(ir) = sqrt(2/pi/Rec(ir).rg/C(1).kr)*exp(-1i*C(1).kr*Rec(ir).rg-pi/4);
        J0L(ir) = sqrt(2/pi/Rec(ir).rg/C(1).kl)*exp(-1i*C(1).kl*Rec(ir).rg-pi/4);
            % sample eigenfunctions at source and receiver depth
            
           if Rec(ir).H==0
               uy0(ir)=C(1).uy(end);
               ux0(ir)=C(1).ux(end);
               uz0(ir)=C(1).uz(end);
           else
               uy0(ir) = interp1(C(1).zz(ib),C(1).uy(ib),Rec(ir).H,'linear');
               ux0(ir) = interp1(C(1).zz(ib),C(1).ux(ib),Rec(ir).H,'linear');
               uz0(ir) = interp1(C(1).zz(ib),C(1).uz(ib),Rec(ir).H,'linear');
            end
       end
       % build Green's function
       for is=1:Src(1).ns
           for ir=1:Rec(1).nr
            % Love waves:
            G1 = J0L(ir);            
            G_L(is,ir).Ghat(ifreq)=uy0(ir)*uyi(is)/...
                (8*C(1).cl*C(1).Ul*C(1).Il(1))*G1;
            G = J0R(ir);
            Gr(1,1) = ux0(ir)*uxi(is)*G;
            Gr(1,2) = -1i*uz0(ir)*uxi(is)*G;
            Gr(2,1) = 1i*ux0(ir)*uzi(is)*G;
            Gr(2,2) = uz0(ir)*uzi(is)*G;
            G_R(is,ir).Ghat(1:2,1:2,ifreq) = Gr/(8*C(1).cr*C(1).Ur*C(1).Ir(1));
         
            
            if strcmp(pol,'rtz')~=1
                C=cos(Rec(ir).az*pi/180); 
                S=sin(Rec(ir).az*pi/180); 
                crap=[Gr(1,1) 0 Gr(1,2);0 0 0;Gr(2,1) 0 Gr(2,2)];
                R=[C S 0;-S C 0;0 0 1];
                crap=R'*crap*R;
                G_R(is,ir).Ghat(1:3,1:3,ifreq) = crap/(8*C(1).cr*C(1).Ur*C(1).Ir(1)); % A&R(2002) 7.147
                G_L(is,ir).Ghat(1:2,1:2,ifreq) = [S2 -S*C;-S*C C^2]*uy0(ir)*uyi(is)/...
                                            (8*C(1).cl*C(1).Ul*C(1).Il(1))*G1; % A&R(2002) 7.146
            end
            
           end
       end
    else %% there are several modes to consider
        
     for imode=1:max_mode
         %% loop over source
       for is=1:Src(1).ns
           if Src(is).H==0
               uyi(is)=C(1).uy(end,imode);
               uxi(is)=C(1).ux(end,imode);
               uzi(is)=C(1).uz(end,imode);
           else
               uyi(is) = interp1(C(1).zz(ib),C(1).uy(ib,imode),Src(is).H,'linear');
               uxi(is) = interp1(C(1).zz(ib),C(1).ux(ib,imode),Src(is).H,'linear');
               uzi(is) = interp1(C(1).zz(ib),C(1).uz(ib,imode),Src(is).H,'linear');
           end
       end
       %% loop over receiver
       for ir=1:Rec(1).nr
        phi(ir)=Rec(ir).az*pi/180; % convert degrees to radian
        J0R(ir) = sqrt(2/pi/Rec(ir).rg/C(1).kr(imode))*exp(-1i*C(1).kr(imode)*Rec(ir).rg-pi/4);
        J0L(ir) = sqrt(2/pi/Rec(ir).rg/C(1).kl(imode))*exp(-1i*C(1).kl(imode)*Rec(ir).rg-pi/4);
            % sample eigenfunctions at source and receiver depth
            
           if Rec(ir).H==0
               uy0(ir)=C(1).uy(end,imode);
               ux0(ir)=C(1).ux(end,imode);
               uz0(ir)=C(1).uz(end,imode);
           else
               uy0(ir) = interp1(C(1).zz(ib),C(1).uy(ib,imode),Rec(ir).H,'linear');
               ux0(ir) = interp1(C(1).zz(ib),C(1).ux(ib,imode),Rec(ir).H,'linear');
               uz0(ir) = interp1(C(1).zz(ib),C(1).uz(ib,imode),Rec(ir).H,'linear');
            end
       end
       % build Green's function
       for is=1:Src(1).ns
           for ir=1:Rec(1).nr
            % Love waves:
            G1 = J0L(ir);            
            G_L(is,ir).Ghat(ifreq,imode)=uy0(ir)*uyi(is)/...
                (8*C(1).cl*C(1).Ul*C(1).Il(1))*G1;
            G = J0R(ir);
            Gr(1,1) = ux0(ir)*uxi(is)*G;
            Gr(1,2) = -1i*uz0(ir)*uxi(is)*G;
            Gr(2,1) = 1i*ux0(ir)*uzi(is)*G;
            Gr(2,2) = uz0(ir)*uzi(is)*G;
            G_R(is,ir).Ghat(1:2,1:2,ifreq,imode) = Gr/(8*C(1).cr*C(1).Ur*C(1).Ir(1));
            
            
            if strcmp(pol,'rtz')~=1
                C=cos(Rec(ir).az*pi/180); 
                S=sin(Rec(ir).az*pi/180); 
                crap=[Gr(1,1) 0 Gr(1,2);0 0 0;Gr(2,1) 0 Gr(2,2)];
                R=[C S 0;-S C 0;0 0 1];
                crap=R'*crap*R;
                G_R(is,ir).Ghat(1:3,1:3,ifreq,imode) = crap/(8*C(1).cr*C(1).Ur*C(1).Ir(1)); % A&R(2002) 7.147
                G_L(is,ir).Ghat(1:2,1:2,ifreq,imode) = [S2 -S*C;-S*C C^2]*uy0(ir)*uyi(is)/...
                                            (8*C(1).cl*C(1).Ul*C(1).Il(1))*G1; % A&R(2002) 7.146
            end
           end
       end
       
     end
       
    end
    disp([num2str(ifreq/length(bigC)*100) ' %'])
end



if strcmp(pol,'rtz')~=1
    for is=1:Src(1).ns
    for ir=1:Rec(1).nr     
        C=cos(Rec(ir).az*pi/180); 
        S=sin(Rec(ir).az*pi/180); 
        crap= G_L(is,ir).Ghat;
        G_L(is,ir).Ghat(1,:) = -S * crap;
        G_L(is,ir).Ghat(2,:) =  C * crap;
        crap= G_R(is,ir).Ghat;
        R=[C S 0;-S C 0;0 0 1];
        for ifreq=1:length(bigC)
           crap(1:2,1:2)= squeeze(G_R(is,ir).Ghat(1:2,1:2,ifreq));
           G_R(is,ir).Ghat(1:2,1:2,ifreq) = R'*crap*R;
        end
    end
    end
                
end
            

if strcmp(typ,'t')==1
for is=1:Src(1).ns
   for ir=1:Rec(1).nr  
        if max_mode==1     
       uu(1:length(bigC))=G_L(is,ir).Ghat(1:length(bigC));
       uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
       G_L(is,ir).G(1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
       for ii=1:2
           for jj=1:2
               uu(1:length(bigC))=G_R(is,ir).Ghat(ii,jj,1:length(bigC));
               uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
               G_R(is,ir).G(ii,jj,1:FT.nwin)=ifft(uu,'symmetric')/FT.dt*fac;
           end
       end
       
        if strcmp(pol,'rtz')~=1     % rotate back from RTD to NED
            C=cos(Rec(ir).az*pi/180); 
            S=sin(Rec(ir).az*pi/180); 
            crap= G_L(is,ir).G;
            G_L(is,ir).Ghat(1,:) = -S * crap;
            G_L(is,ir).Ghat(2,:) =  C * crap;
            crap= G_R(is,ir).Ghat;
            R=[C S 0;-S C 0;0 0 1];
            for ifreq=1:FT.nwin
               crap(1:2,1:2)= squeeze(G_R(is,ir).G(1:2,1:2,ifreq));
               G_R(is,ir).G(1:3,1:3,ifreq) = R'*crap*R;
            end
            
            
        end
        else
            for imode=1:max_mode
               uu(1:length(bigC))=G_L(is,ir).Ghat(1:length(bigC));
               uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
               G_L(is,ir).G(1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt*fac;
               for ii=1:2
                for jj=1:2
                   uu(1:length(bigC))=G_R(is,ir).Ghat(ii,jj,1:length(bigC));
                   uu(length(bigC)+1:FT.nwin)=0;uu(1)=0;
                   G_R(is,ir).G(ii,jj,1:FT.nwin,imode)=ifft(uu,'symmetric')/FT.dt*fac;
                end
               end
            end       
            
            
            
        end
   end
end
end


return