%% AIR  - SOLID with 2 solid layers over half space
% Example of a simple gradient over half space. 

% Marine Denolle (09/2014)
addpath '../../sub'
clear all;
close all;clc
mkdir('plots') 
mkdir('mat')


%% TIME /FREQUENCY:
typ = 't' ;         % 't': enters time series length (win= in s) and rate (rate=1/s)
                    % 'f': enters frequency series length (win=in Hz) and rate
                    % (1/Hz)
win = 500;          % (in s) window length
rate = 0.05;       % (in 1/s) sampling rate (10Hz)
FT = make_ft(win,rate,typ);
FT.Fmax=2;         % calculate until Fmax if < Nyquist
max_mode=1;         % maximum number of surface-wave modes to ouput
freq=1./(1:10);     % plots eigenfunctions for periods 1 -> 10s


%% MEDIUM
H = [0 5 20] ;     % (km) vector of interface depths 
beta = 3*[0.5 1];   % (km/s) S-wavespeed 
alpha = 5.4*[0.5 1];%sqrt(3)*beta; % (km/s) P-wavespeed 
rho = 2.27*[0.5 1];    % (kg/dm^3) density
Dmax=30*max(beta)/FT.df;      % (km) "half space" effective thickness, matters for long period
typ = {'lin','cst'};          % constant within layers: 'cst', all gradients 'grad'
MED = make_layers(H,Dmax,alpha,beta,rho,typ);
!mv *.ps plots


%% RESOLUTION:
% NR : resolution structure for adaptive sampling (see Section 2.3)
% typ = 'lin';       % bounded linear(ramp) adaptive sampling 
%     lambda=[2 10]; % H/lambda domains of sampling laws
%     a1=[0.5 1.5];  % ramp lower corner
%     a2=[2 8];      % ramp upper corner
% 	N1=[20 150];    % lower bound # points
% 	N2 = [40 100];  % upper bound # points
%     NR = make_npt_colloc_lin(lambda,N1,N2,a1,a2);
   

NR.typ='rec';       % choose recursive
NR.N1 = 30 ;           % lower bound for number of points / layer
NR.N2 = 100 ;        % maximum number of points / layers


%% GESC
 bigC = gesc(MED,FT,NR,max_mode,freq); % plot eigenfunctions at specific
% frequencies
 save(['./mat/air2solid.mat'],'bigC','NR','MED','FT','-v7.3')
%load mat/air2solid
%% get depth correction
for HH=1:10 % loop over potential source depths from 1 to 10 km
for i=1:length(bigC)
    C=bigC(i).a;
    if isempty(C);continue;end
    [~,ib]=unique(C(1).zz);
      FREQ(i)=C(1).omega/(2*pi);
      R10 = interp1(C(1).zz(ib),C(1).ux(ib),0,'linear');
      R20 = interp1(C(1).zz(ib),C(1).uz(ib),0,'linear');
      L10 = interp1(C(1).zz(ib),C(1).uy(ib),0,'linear');
      R1(i,HH) = interp1(C(1).zz(ib),C(1).ux(ib),HH,'linear')/R10;
      R2(i,HH) = interp1(C(1).zz(ib),C(1).uz(ib),HH,'linear')/R20;
      L1(i,HH) = interp1(C(1).zz(ib),C(1).uy(ib),HH,'linear')/L10;
end
end
ik=find(FREQ~=0);
figure
subplot(311)
plot(FREQ(ik),R1(ik,:));grid on
set(gca,'fontsize',14);xlabel('Freq (Hz)');title('R1 depth correction');
subplot(312)
plot(FREQ(ik),R2(ik,:));grid on
set(gca,'fontsize',14);xlabel('Freq (Hz)');title('R2 depth correction');
subplot(313)
plot(FREQ(ik),L1(ik,:));grid on
set(gca,'fontsize',14);xlabel('Freq (Hz)');title('L1 depth correction');


% have fun with the rest of the code.

%% Get dispersion curves
[freq,Cr,Cl,Ur,Ul]=get_disp(bigC,max_mode);
ik=find(Cr(:,1)~=0);
figure(5)
for i=1:max_mode
    subplot(211)
    plot(freq(ik),Cr(ik,i),'r',freq(ik),Cl(ik,i),'b');
    grid on;hold on;set(gca,'Fontsize',14);xlabel('Frequency (Hz)');
    legend('Phase R','Phase L')
    subplot(212)
    plot(freq(ik),Ur(ik,i),'r',freq(ik),Ul(ik,i),'b');
    grid on;hold on;set(gca,'Fontsize',14);xlabel('Frequency (Hz)');
    legend('Group R','Group L')
end
xlim([0 1])
print('-dpsc',['./plots/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ) '_disp.ps']);
% save(['./mat/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ)],'bigC','NR','MED','FT',...
%     'Cr','Cl','Ur','Ul','-v7.3');
% 
%% Get surface Rayleigh-wave particle motion
Z=(0:10);               % depth at which calculating particle motion
[freq,ux,uz]=get_sw_pm(bigC,max_mode,Z);
figure(6)
for iz=1:1;%length(Z)
for i=1:max_mode
    plot(freq,ux(iz,:,i)./uz(iz,:,i));
    xlim([0 2]);ylim([-1 1])
    grid on;hold on;set(gca,'Fontsize',14);xlabel('Frequency (Hz)');
    title('Ellipticity')
end
end
print('-dpsc',['./plots/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ) '_sw_pm.ps']);
% save(['./mat/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ)],'bigC','NR','MED','FT',...
%    'Cr','Cl','Ur','Ul','ux','uz','-v7.3');


%% Get surface H/V ratio
[freq,HV]=get_hv(bigC,max_mode);
figure(6)
for i=1:max_mode
    plot(freq,HV(:,i));
    grid on;hold on;set(gca,'Fontsize',14);xlabel('Frequency (Hz)');
    title('H/V Ration')
end
print('-dpsc',['./plots/MED_' int2str(MED(1).Nl-1) 'l_' char(MED(1).typ) '_sw_pm.ps']);


%% Get SW Green's functions
% define sources:
Src(1).ns=1; % source different sources
Src(1).H =0;% Src(2).H=6; % two different source depths
% choose receivers:
Rec(1).nr=100;
Rgmax=400;
Rgmin=50;
N1=20;
N2=20;
for i=1:N1%Rec(1).nr
    for ii=1:N2
        jj=N1*(i-1)+ii;
    Rec(jj).az = 360*rand(1); 
    Rec(jj).rg = (Rgmax-Rgmin)/(N2)*i+Rgmin;
    Rec(jj).H=0;
    end
end
Rec(1).nr=N1*N2;
[G_L,G_R] = get_sw_gf_farfield(FT,bigC,Src,Rec,max_mode,'t','rtz');
for is=1:Src(1).ns
%      figure
   for ir=1:Rec(1).nr      
%        sina=sin(Rec(ir).az*pi/180);sinb=sina;
%        cosa=cos(Rec(ir).az*pi/180);cosb=cosa;
%      for it=1:FT.nwin
%          crap=G_L(is,ir).G(1:2,1:2,it) ;
%         TT(ir,it) = sina*sinb*crap(1,1)-sina*cosb*crap(1,2)-cosa*sinb*crap(2,1)+cosa*cosb*crap(2,2);
%          crap=G_R(is,ir).G(1:3,1:3,it) ;
%         RR(ir,it) = cosa*cosb*crap(1,1)+cosa*sinb*crap(1,2)+sina*cosb*crap(2,1)+sina*sinb*crap(2,2);
%         RZ(ir,it) = cosa*crap(1,3)+sina*crap(2,3);
%         ZR(ir,it) = -sina*crap(1,3)+cosa*crap(2,3);
%         disp([RZ(ir,it) ZR(ir,it)])
%         ZZ(ir,it) = crap(3,3);
%      end
     TT(ir,:)=G_L(is,ir).G;
     RR(ir,1:FT.nwin)=squeeze(G_R(is,ir).G(1,1,:));
     RZ(ir,1:FT.nwin)=squeeze(G_R(is,ir).G(1,2,:));
     ZR(ir,1:FT.nwin)=squeeze(G_R(is,ir).G(2,1,:));
     ZZ(ir,1:FT.nwin)=squeeze(G_R(is,ir).G(2,2,:));
 
%      subplot(3,3,5)
%      plot(FT.tt(1:FT.nwin)-win/2,TT(ir,:)*1E16+Rec(ir).rg,'b');hold on;grid on
%      subplot(3,3,1)
%      plot(FT.tt(1:FT.nwin)-win/2,RR(ir,:)*1E16+Rec(ir).rg,'b');hold on;grid on
%      subplot(3,3,3)
%      plot(FT.tt(1:FT.nwin)-win/2,RZ(ir,:)*1E16+Rec(ir).rg,'b');hold on;grid on
%      subplot(3,3,7)
%      plot(FT.tt(1:FT.nwin)-win/2,ZR(ir,:)*1E16+Rec(ir).rg,'b');hold on;grid on
%      subplot(3,3,9)
%      plot(FT.tt(1:FT.nwin)-win/2,ZZ(ir,:)*1E16+Rec(ir).rg,'b');hold on;grid on
%      
%      
%         
%      pause(1)
   end
end

 save('./mat/air2solid_GF_az.mat','ZZ','TT','ZR','RZ','RR','Rec','FT','G_L','G_R')
 
%% Get surface-wave wavefield in 2D
% define sources:
Src(1).ns=1; % source different sources
Src(1).H = 6;% Src(2).H=6; % two different source depths
% Source is going to follow a Brune spectrum:
Src(1).Mw=6;
Ds=3E6 ; % average stress drop 3MPa
Src(1).Mo=10^(3/2*Src(1).Mw+9); % (Nm)
fc = 1500*(Ds/norm(8.5*Src(1).Mo,2))^(1/3); % corner frequency
tau = 1/2*1/fc; % pulse width 1/(2*fc)
Src(1).Shat(2:length(FT.omega)) = Src(1).Mo.*exp(-1i*FT.omega(2:end).*3/2*tau)./(1+(FT.omega(2:end)/(2*pi*fc)).^2)./(1i*FT.omega(2:end));
Src(1).Shat(length(FT.omega)+1:FT.nwin)=0;
Src(1).Stf = ifft(Src(1).Shat,'symmetric')/FT.dt; % conjugate due to sign convention on FFT

Src(1).M=[0 1 0;1 0 0;0 0 0];
Src(1).H = 6;% Src(2).H=6; % two different source depths
% choose receivers in array of 50 by 30 receivers on a horizontal plane
% Nx=5;
% Rmax=300; %km
% R=linspace(0,Rmax,Nx);
% H=linspace(0,40,Nx);
% Rec(1).nr = length(R)*length(H);
% for ir1=1:length(R)
%     for ir2=1:length(H)
%         Rec(Nx*(ir1-1)+ir2).H=H(ir2);
%         Rec(Nx*(ir1-1)+ir2).rg=R(ir1);
%         Rec(Nx*(ir1-1)+ir2).az=0;
%     end
% end
% [G_L,G_R] = get_sw_gf(FT,bigC,Src,Rec,max_mode,'t');
% [nerr]=plot_sw_gf_2D(FT,G_L,G_R,Src,Rec,max_mode,'xz','rtz');
 

%% To Get SW response to single point force source, multiply Green tensor with single force


%% Get SW response to double couple
[Ldc,Rdc] = get_sw_dc_farfield(FT,bigC,Src,Rec,max_mode,'t','rtz');
figure
[a,b]=butter(2,[1/20 2]/10);
for ir=1:Rec(1).nr      
 T(ir,:)=filtfilt(a,b,Ldc(is,ir).U);
 R(ir,1:FT.nwin)=filtfilt(a,b,squeeze(Rdc(is,ir).U(1,:)));
 Z(ir,1:FT.nwin)=filtfilt(a,b,squeeze(Rdc(is,ir).U(2,:)));

%  subplot(1,3,1)
%  plot(FT.tt(1:FT.nwin),T(ir,:)*1000+Rec(ir).rg,'b');hold on;grid on
%  subplot(1,3,2)
%  plot(FT.tt(1:FT.nwin),R(ir,:)*1000+Rec(ir).rg,'b');hold on;grid on
%  subplot(1,3,3)
%  plot(FT.tt(1:FT.nwin),Z(ir,:)*1000+Rec(ir).rg,'b');hold on;grid on
% 
%  
% 
%  pause(1)
end
save('./mat/air2solid_GF_az.mat','Src','ZZ','TT','ZR','RZ','RR','Rec','FT','G_L','G_R','R','T','Z')
 

%% get dynamic strain and stress at any depth
[EE,SS] = get_strain_stress(MED,FT,bigC,R,T,Z,5,'rtz');
