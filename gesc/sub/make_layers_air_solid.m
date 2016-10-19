%% Create Medium structure MED using layers
function MED = make_layers_air_solid(H,Dmax,alpha,beta,rho,typ)

%% input
% H: (in km) vector of layers depths starting with depth = 0 km (at the
% surface) and positive downward
% Dmax: input thickness of half space: bottom layer may be taken thick to
% simulate a half space: take it large for broadband caluclations
% alpha: (in km/s) vector of P-wave speed (Vp) of length *length(H)-1*
% beta: (in km/s) vector of S-wave speed (Vs) of length *length(H)-1*
% rho: (in kg/dm3) vector of density of length *length(H)-1*
% typ: vector of {'cst'} (homogeneous layers), {'lin'} (constant gradient layers)

%% output
% MED: structure with medium properties. GESC requires organizing medium
% from bottom to top, and for elastic properties *within* layers to match
% properly boundary conditions

% the GESC code asks for the elastic properties from the bottom of the
% MEDium to the top:
% MED is a vector (length = Nl number of layers) of structures
% MED(1) contains depth, Vs, Vp, Rho and depth of interfaces for entire
% MEDium
% MED(1) corresponds to the half space
% MED(2:Nl) contains depth, Vs, Vp, Rho *within* the layers, from bottom to
% top

% Note: to choose a Poisson medium, poiss=1 and alpha = sqrt(3)beta
%       to choose homogeneous layers, use typ = 'cst'
%       to choose constant gradient layers, use typ = 'lin'
% Marine Denolle (04/2014)


% check if medium is Poisson (i.e. if Vp = sqrt(3)*Vs)
MED(1).poiss = 0; % not a Poisson medium
if all(alpha==sqrt(3)*beta) ; MED(1).poiss = 1;end

% create medium in vectors
z=[];a=[];b=[];r=[];
alpha = [alpha alpha(end)];
beta = [beta beta(end)];
rho = [rho rho(end)];
for i=2:length(H)
    disp(typ{i-1})
    if strcmp(typ{i-1},'cst')==1
        z = [linspace(H(i),H(i-1),50) z]; % depth vector
        a = [alpha(i-1)*ones(1,50) a]; % alpha / P wave velocity vector
        b = [beta(i-1)*ones(1,50) b]; % beta / S wave velocity vector
        r = [rho(i-1)*ones(1,50) r]; % rtho /  density vector
    elseif strcmp(typ{i-1},'lin')==1
        z = [linspace(H(i),H(i-1),50) z]; % depth vector
        a = [linspace(alpha(i),alpha(i-1),50) a]; % alpha / P wave velocity vector
        b = [linspace(beta(i),beta(i-1),50) b]; % beta / S wave velocity vector
        r = [linspace(rho(i),rho(i-1),50) r]; % rho /  density vector
    else
        disp('Enter type of layer poperly')
        return
    end
end


% remove overlapping point (2 points have the same depth H1, H1)
[ia,ib]= unique(z);
% outpout data in MED
MED(1).z=ia ; MED(1).z = z(ib);
MED(1).alpha = a(ib);
MED(1).beta = b(ib);
MED(1).rho = r(ib);
MED(1).inter = H;
MED(1).Nl=length(H);
% break into layers: elastic properties have to be contained within layer
for i=2:MED(1).Nl
    ik=find(MED(1).inter(i-1)<MED(1).z&MED(1).inter(i)>MED(1).z);
    MED(MED(1).Nl-i+2).betal = flipud(MED(1).beta(ik));
    MED(MED(1).Nl-i+2).alphal = flipud(MED(1).alpha(ik));
    MED(MED(1).Nl-i+2).rhol = flipud(MED(1).rho(ik));
    MED(MED(1).Nl-i+2).zz = flipud(MED(1).z(ik));
    MED(MED(1).Nl-i+2).zz(1) = MED(1).inter(i);
    MED(MED(1).Nl-i+2).zz(end) = MED(1).inter(i-1);
    MED(MED(1).Nl-i+2).typ=typ{i-1};
end
MED(1).inter = fliplr(MED(1).inter);
MED(1).Dmax=Dmax;

% output velocity/density profiles
figure(102)
plot(b,z,'b-o',a,z,'r-o',r,z,'k-o','Linewidth',2);set(gca,'Fontsize',14);grid on;axis ij
ylim([0 max(H)]);ylabel('Depth (km)');legend('Vs (km/s)','Vp (km/s)','density (kg/dm^3)')
title(['Depth profile ' int2str(MED(1).Nl) ' layers']);
set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 8 7],'PaperPositionMode','manual');
print('-dpsc',['./eig_' int2str(MED(1).Nl-1) 'l.profile.ps'])

return
