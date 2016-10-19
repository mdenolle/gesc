% Function to extract Rayleigh-wave particle motion at depth Z
function [freq,ux,uz]=get_sw_pm(bigC,max_mode,Z)
%% INPUTS:
% bigC: gesc eigen structure
% max_mode: maximum modes to evaluate the particle motion on
% Z : vector of depths (in km) to evaluate the particle motion on
%% OUTOUTS:
% ux: horizontal Rayleigh waves at all frequencies, modes and depth
% requested
% uz: vertical Rayleigh waves at all frequencies, modes, and depth
% requested

for i=1:length(bigC)
    C=bigC(i).a;
    if isempty(C);continue;end
    freq(i)=C(1).omega/(2*pi);
    [~,ib]=unique(C(1).zz);
    for iz=1:length(Z)
        for imode = 1:max_mode
         ux(iz,i,imode) = interp1(C(1).zz(ib),C(1).ux(ib,imode),Z(iz),'spline');
         uz(iz,i,imode) = interp1(C(1).zz(ib),C(1).uz(ib,imode),Z(iz),'spline');
        end
    end
end
ux=squeeze(ux);uz=squeeze(uz); % if only one depth requested, reduce matrix of dimension 3 to 2.
return

