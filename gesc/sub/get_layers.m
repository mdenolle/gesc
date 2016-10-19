% function to pick interfaces from velocity profiles
% INPUT:
% z: (km) depth vector 
% a: (km/s) P-wave speed
% b: (km/s) S-wave speed
% r: (kg/dm3) density
% OUTPUT:
% HH: (km) vector of interfaces, descending order, will end with the
% surface at zero depth automatically.
% USAGE:
%  ------- zoom in ----------------------
%    2 left clicks on left mouse button
%  ------- zoom out ---------------------
%    press "z" on keyboard (channel 112)  
%  ------- pick interface ---------------
%    press "p" on keyboard (channel 112)   
%    left click on the according depth (only depth matter)
%  ------- delete previous pick ---------
%    press "d" on keyboard (channel 100)
%  ------- delete all picks -------------
%    press "e" on keyboard (channel 101)
%  ------- done with picks -------------
%    press "q" on keyvboard (channel 113)
% Marine Denolle (04/08/2014)
function HH = get_layers(z,a,b,r)

close all

HH=[];
figure(10000)
plot(b,z,'b-o',a,z,'r-o',r,z,'k-o','Linewidth',2);set(gca,'Fontsize',14);grid on;axis ij
hold on
ylim([0 max(z)]);ylabel('Depth (km)');legend('Vs (km/s)','Vp (km/s)','density (kg/dm^3)')

NL=0; % number of interfaces;

b=0; % keys pressed or mouse clicks
while (b~=113)

	[~,z1,b]=ginput(1);	% upper layer
    same=true;
    
    while (same)
        switch b
            case 1 % left click
                [~,z2,b]=ginput(1);	% lower layer
                if (b==1) && z2>z1
                     ylim([z1 z2])
                end
            case 122 % right click , zoom out
                ylim('auto');xlim('auto')
            case 112
                NL=NL+1;
                [~,z2,b]=ginput(1);	% lower layer
                HH(NL) = z2;
                h{NL}=plot(linspace(0,10,10),z2*ones(10,1),'r');
            case 100 % delete the previously picked interface
                disp(NL)
                delete(h{NL})
                HH(NL)=0;
                NL=NL-1;
            case 101 % delete all picks
                for i=1:NL
                    disp('deleting interfaces')
                    delete(h{NL})
                end
                NL=0;
                HL=[];
            case 113 % quit
                same=false;
        end
	if (same),[~,z1,b]=ginput(1);end	%clic ou frappe suivante
    end
end
close all
HH=[sort(HH,'descend') 0];
return