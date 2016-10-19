%% Create Medium structure MED using layers
function MED = read_medium(fname,Hmax,Dmax)
% reads in fname the velocity and density profiles
%% input:
%  fname: file name (advise absolute path)
%  Hmax: (km) maximum depth to read from fname
%  Dmax: (km) depth of the "half space"


% the file should have the following structure:
% z_1 Vp_1 Vs_1 Rho_1
%  .   .    .     .
%  .   .    .     .
%  .   .    .     .
% GESC requires the wavespeed to be (km/s), density (kg/dm3) and depth (km)
% and will verify so
% note: the depths have to be in km, the scripts verifies this based on the
% bottom depth: if > mean(depth)>1E4 (in average more than 10,000) then convert m to km. 
%% ouput:
% MED: GESC structure for medium


%% get short name for fname:
toto=char(fname);
for i=length(toto):-1:1
    if strcmp(toto(i:i),'/')==1
        toto2 = toto(i+1:end);
        break;
    end
end
toto=char(toto2);
for i=length(toto):-1:1
    if strcmp(toto(i:i),'.')==1
        truename = toto(1:i-1);
        break;
    end
end
disp(truename)
MED(1).name = truename;
MED(1).Dmax=Dmax;


%% read in data:
x=importdata(fname);

% convert to km,km/s and to kg/fm3
fracv=1;fracr=1;fracz=1;
if mean(x(1,:))>1E4;fracz=1E-3;end
if x(1,2)>50;fracv=1E-3;end
if x(1,end)>50;fracr=1E-3;end
if x(end,1) <= Hmax 
    i_end=length(x(:,1));
else
    [~,i_end]=min(abs(x(:,1)-Hmax));
end
% format depth formats
if x(1,1)>x(end,1) %decreasing depths
    MED(1).z = flipdim(x(1:i_end,1),1)*fracz;
    MED(1).alpha = flipdim(x(1:i_end,2),1)*fracv;
    MED(1).beta = flipdim(x(1:i_end,3),1)*fracv;
    MED(1).rho = flipdim(x(1:i_end,4),1)*fracr;
else % increasing depth
    MED(1).z = x(1:i_end,3)*fracz;
    MED(1).alpha = x(1:i_end,2)*fracv;
    MED(1).beta = x(1:i_end,3)*fracv;
    MED(1).rho = x(1:i_end,4)*fracr;
end

%% plots
plot(MED(1).beta,MED(1).z,'b-o',MED(1).alpha,MED(1).z,'r-o',...
    MED(1).rho,MED(1).z,'k-o','Linewidth',2);set(gca,'Fontsize',14);grid on;axis ij
ylim([0 max(MED(1).z)]);ylabel('Depth (km)');legend('Vs (km/s)','Vp (km/s)','density (kg/dm^3)')
title(['Depth profile from '  char(truename)]);
% set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 8 7],'PaperPositionMode','manual');
% print('-dpsc',['./PLOTS/eig_' char(truename) '_profile.ps'])

disp('Please check where you need to impose layers')
disp('Remember that:')
disp('     - interfaces may not be seen as such for given frequencies')
disp('     - too many layers will scale UP the system')
disp('     - smoothly varying properties should stay within one layer')
disp('You are about to pick the layers using get_layers.m interactive subroutine')
disp('Look carefully now the layers ... and click enter when done')
pause
close all


HH=get_layers(MED(1).z,MED(1).alpha,MED(1).beta,MED(1).rho);
dh=abs(MED(1).z(2)-MED(1).z(1));
% clean up
for i=1:length(HH)
    if min(abs(HH(i)-MED(1).z))<dh/4 % if chosen point is less than 1/2 spacing, then pick point
        [~,ikk]=min(abs(HH(i)-MED(1).z));
        HH(i)=MED(1).z(ikk);
    end
end

%% layering:
MED(1).inter = [HH Hmax]; 
MED(1).Nl=length(MED(1).inter);
for i=2:MED(1).Nl
    ik=find(MED(1).inter(i-1)<=MED(1).z&MED(1).inter(i)>=MED(1).z);
    MED(MED(1).Nl-i+2).betal = flipud(MED(1).beta(ik));
    MED(MED(1).Nl-i+2).alphal = flipud(MED(1).alpha(ik));
    MED(MED(1).Nl-i+2).rhol = flipud(MED(1).rho(ik));
    MED(MED(1).Nl-i+2).zz = flipud(MED(1).z(ik));
    MED(MED(1).Nl-i+2).zz(1) = MED(1).inter(i);
    MED(MED(1).Nl-i+2).zz(end) = MED(1).inter(i-1);
end
MED(1).inter = fliplr(MED(1).inter);
save(['./MAT/MED_' char(truename) '.mat'],'MED')

plot(MED(1).beta,MED(1).z,'b-o',MED(1).alpha,MED(1).z,'r-o',...
    MED(1).rho,MED(1).z,'k-o','Linewidth',2);set(gca,'Fontsize',14);grid on;axis ij
hold on
for i=1:length(HH)
    plot(linspace(0,10,10),HH(i)*ones(10,1),'r');
end
xlim([0 1+max(MED(1).alpha)])
ylim([0 max(MED(1).z)]);ylabel('Depth (km)');legend('Vs (km/s)','Vp (km/s)','density (kg/dm^3)')
title(['Depth profile from '  char(truename)]);
set(gcf,'PaperUnits','inches','PaperPosition',[0.25 0.25 8 7],'PaperPositionMode','manual');
print('-dpsc',['./PLOTS/eig_' char(truename) '_profile.ps'])


return