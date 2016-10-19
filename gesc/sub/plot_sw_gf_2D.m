% plot the 2D wavefield of the Green tensor components
function fname = plot_sw_gf_2D(FT,G_L,G_R,Src,Rec,max_mode,typ,coo)
% typ = 'xy' for horizontal
% coo: 'rtz' or 'zrt'

fname='';
if strcmp(coo,'rtz')==1
    chan={'RR','RT','RZ','TR','TT','TZ','ZR','ZT','ZZ'};
elseif strcmp(coo,'nez')==1
    chan={'NN','NE','NZ','EN','EE','EZ','ZN','ZE','ZZ'};
else
    disp('Impose correct coordinate system: nez or rtz?')
    return
end
if length(G_L(1,:))~=Rec(1).nr||length(G_R(1,:))~=Rec(1).nr||...
   length(G_L(:,1))~=Src(1).ns||length(G_R(:,1))~=Src(1).ns;
disp('wrong sizes');return;end


cc=[cmap([0 0 1],50,40,0);flipud(cmap([1 0 0],50,40,0))];
cscale= linspace(-1,1,length(cc(:,1)));

fname=['./plots/animations/GF_2D/' char(typ) '/Nsource_' int2str(1) ...
    '_Nreceiver_' int2str(Rec(1).nr)];
mkdir(fname)
if strcmp(typ,'xy')==1
    % the plot will be projected on a horizontal plane
    % get mesh
    rg=zeros(Rec(1).nr,1);az=rg;
    for i=1:Rec(1).nr
        rg(i)=Rec(i).rg;
        az(i)=Rec(i).az;
    end
    ii=find(rg==min(rg));Nx=length(ii);GG=zeros(Nx);
    for imode=1:max_mode

    for is=1:Src(1).ns
      fdir = [char(fname) '/SOURCE_' int2str(is)];mkdir(fdir)
        if strcmp(coo,'nez')==1 % sum all components in the NEZ framework
            for ii=1:3;for jj=1:3;
            for it=1:2:FT.nwin
            time = FT.tt(it);
            for ir=1:Rec(1).nr
             if mod(ir,Nx)==0;
                i1=floor(ir/Nx);i2=Nx;
             else
                 i1=floor(ir/Nx)+1;i2=mod(ir,Nx);
             end
               if ii <=2 && jj <=2
              GG(i1,i2) = G_L(is,ir).G(ii,jj,it)+G_R(is,ir).G(ii,jj,it);
                else
              GG(i1,i2) = G_R(is,ir).G(ii,jj,it);
                end
              X(i1)=Rec(ir).rg;
              Y(icrapi2)=Rec(ir).az*pi/180;
            end
            [YY,XX]=meshgrid(Y,X);
            [AA,BB]=pol2cart(YY,XX);
            pcolor(AA,BB,GG);
            if time<=2
                max1=max(max(abs(GG)));
            end
            colormap(cc);shading flat;
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,0,'rp','Markersize',10,'Linewidth',3);axis square
            ylabel('Distance (km)');title(['channel ' char(chan(3*(ii-1)+jj)) ' ' num2str(FT.tt(it)) ' seconds'])
            caxis(max1*[-1 1]);
            if it==1;pause(0.25);end
            [char(fdir) '/' char(chan(3*(ii-1)+jj)) '_' int2str(it) '.jpg']
            print('-djpeg','-r150',[char(fdir) '/' char(chan(3*(ii-1)+jj)) '_' int2str(it) '.jpg'])
            end % end loop over time
            end;end % end loops over channels
        else   % rotate into RTZ !!!!
             
            for it=1:2:FT.nwin
            time = FT.tt(it);
            crap=zeros(3);
            for ir=1:Rec(1).nr
                
             for ii=1:3;for jj=1:3;
              if ii <=2 && jj <=2
               crap(ii,jj) = G_L(is,ir).G(ii,jj,it)+G_R(is,ir).G(ii,jj,it);
              else
               crap(ii,jj) = G_R(is,ir).G(ii,jj,it);
              end
              end;end
              crap1=rot_nez2rtz(crap,Rec(ir).az);
              GG(floor(ir/Nx)+1,mod(ir,Nx)+1,1:9,it)=reshape(crap1,9,1);     
              X(floor(ir/Nx)+1)=Rec(ir).rg;
              Y(mod(ir,Nx)+1)=Rec(ir).az*pi/180;
            end
            end
            [YY,XX]=meshgrid(Y,X);
            [AA,BB]=pol2cart(YY,XX);
            crap=zeros(Nx);
            for jj=1:9
              for it=1:4:FT.nwin
               crap(1:Nx,1:Nx)=GG(1:Nx,1:Nx,jj,it);
               pcolor(AA,BB,crap)
               if time<=2
                max1=max(max(abs(crap)));
               end
             colormap(cc);shading interp;
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,0,'rp','Markersize',10,'Linewidth',3);axis square
            ylabel('Distance (km)');title(['channel ' char(chan(3*(ii-1)+jj)) ' ' num2str(FT.tt(it)) ' seconds'])
            caxis(max1*[-1 1]);
            if it==1;pause(0.25);end
            [char(fdir) '/' char(chan(3*(ii-1)+jj)) '_' int2str(it) '.jpg']
            print('-djpeg','-r150',[char(fdir) '/' char(chan(jj)) '_' int2str(it) '.jpg'])
               end % end loop over time
            end % end loops over channels
            
        end
            
            
            
    end
    end
else
      % the plot will be projected on a horizontal plane
    % get mesh
    for i=1:Rec(1).nr
        rg(i)=Rec(i).rg;
        HH(i)=Rec(i).H;
    end
    ii=find(rg==min(rg));Nx=length(ii);
    ii=find(HH==min(HH));Ny=length(ii);GG=zeros(Nx,Ny);
    for imode=1:max_mode

    for is=1:Src(1).ns
      fdir = [char(fname) '/SOURCE_' int2str(1)];mkdir(fdir)
        % arrange such that N = R, E=T and Z=Z
            for ii=1:3;for jj=1:3;
            for it=1:4:FT.nwin
            time = FT.tt(it);
             for ir=1:Rec(1).nr
             if mod(ir,Nx)==0;
                i1=floor(ir/Nx);i2=Nx;
             else
                 i1=floor(ir/Nx)+1;i2=mod(ir,Nx);
             end
        
              if ii <=2 && jj <=2
              GG(i1,i2) = G_L(is,ir).G(ii,jj,it)+G_R(is,ir).G(ii,jj,it);
                else
              GG(i1,i2) = G_R(is,ir).G(ii,jj,it);
                end
              X(i1)=Rec(ir).rg;
              Y(i2)=Rec(ir).H;
             end
            [YY,XX]=meshgrid(Y,X);
            pcolor(XX,YY,GG);axis ij
            if time <=10
                max1=max(max(abs(GG)));
            end
            colormap(cc);shading interp;axis equal;axis image
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,0,'rp','Markersize',10,'Linewidth',3);
            ylabel('Depth (km)');title(['channel ' char(chan(3*(ii-1)+jj)) ' ' num2str(FT.tt(it)) ' seconds'])
            caxis(max1/10*[-1 1])
            if it==1;pause(0.25);end
            [char(fdir) '/R_Z_' char(chan(3*(ii-1)+jj)) '_' int2str(it) '.jpg']
            print('-djpeg','-r150',[char(fdir) '/' char(chan(3*(ii-1)+jj)) '_' int2str(it) '.jpg'])
            end % end loop over time
            end;end % end loops over channels
      
            
        end
    end
    
    
    
end
           