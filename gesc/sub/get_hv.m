% subroutine to get the  Horizontal/Vertical ratio.
% Note that the horizontal components include both Rayleigh and Love waves
function [freq,HV]=get_hv(bigC,max_mode)
%% INPUTS:
% freq : frequency vector (in Hz) to estimate the HV ratio on
% bigC: gesc eigen structure
% max_mode: maximum modes to evaluate the HV on
%% OUTOUTS:
% HV: HV ratio evaluated at all frequencies and all modes requested.


for i=1:length(bigC)
    C=bigC(i).a;
    if isempty(C);continue;end
    freq(i)=C(1).omega/(2*pi);
    for imode=1:max_mode
      HV(i,imode) = abs(C(1).ux(end,imode)/C(1).uz(end,imode));
    end
end
return
