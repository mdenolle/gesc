% subroutine to get the dispersion curves
function [freq,cr,cl,Ur,Ul]=get_disp(bigC,max_mode)


for i=1:length(bigC)
    C=bigC(i).a;
    if isempty(C);continue;end
    freq(i)=C(1).omega/(2*pi);
   for imode = 1:length(C(1).cr)
      cr(i,imode) = C(1).cr(imode);
      cl(i,imode) = C(1).cl(imode);
      Ur(i,imode) = C(1).Ur(imode);
      Ul(i,imode) = C(1).Ul(imode);
end
end


return