% plot the 2D wavefield of the Green tensor components
function fname = plot_sw_2D(fdir,FT,G_L,G_R,Src,Rec,max_mode,typ,coo)
% typ = 'xy' for horizontal
% coo: 'rtz' or 'zrt'

fname='';
if strcmp(coo,'rtz')==1
    chan={'R','T','Z'};
elseif strcmp(coo,'nez')==1
    chan={'N','E','Z'};
else
    disp('Impose correct coordinate system: nez or rtz?')
    return
end
if length(G_L(1,:))~=Rec(1).nr||length(G_R(1,:))~=Rec(1).nr||...
   length(G_L(:,1))~=Src(1).ns||length(G_R(:,1))~=Src(1).ns;
disp('wrong sizes');return;end


cc=[cmap([0 0 1],50,40,0);flipud(cmap([1 0 0],50,40,0))];
cscale= linspace(-1,1,length(cc(:,1)));

fname=[char(fdir) '/plots/animations/GF_2D/' char(typ) '/Nsource_' int2str(Src(1).ns) ...
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
    X=zeros(Nx,1);Y=X;
    for imode=1:max_mode

    for is=1:Src(1).ns
      fdir = [char(fname) '/SOURCE_' int2str(is)];mkdir(fdir)
      if strcmp(coo,'nez')==1 % sum all components in the NEZ framework
        for ii=1:3
          for it=1:2:FT.nwin
            for ir=2:Rec(1).nr
             if mod(ir,Nx)==0;
                i1=floor(ir/Nx);i2=Nx;
             else
                 i1=floor(ir/Nx)+1;i2=mod(ir,Nx);
             end
                if ii <=2
              GG(i1,i2) = G_L(is,ir).U(ii,it,imode)+G_R(is,ir).U(ii,it,imode);
                else
              GG(i1,i2) = G_R(is,ir).U(ii,it,imode);
                end
              X(i1)=Rec(ir).rg;
              Y(i2)=Rec(ir).az*pi/180;
            end
            [YY,XX]=meshgrid(Y,X);
            [AA,BB]=pol2cart(YY,XX);
            pcolor(AA,BB,GG);
              if it*FT.dt <=5
                max1=max(max(abs(GG)));
              end
            colormap(cc);shading interp;
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,0,'rp','Markersize',10,'Linewidth',3);axis square
            ylabel('Distance (km)');title(['channel ' char(chan(ii)) ' ' num2str(FT.tt(it)) ' seconds'])
            caxis(max1/5*[-1 1])
            if it==1;pause(0.25);end
            [char(fdir) '/' char(chan(ii)) '_' int2str(it) '.jpg']
            print('-djpeg','-r150',[char(fdir) '/' char(chan(ii)) '_' int2str(it) '.jpg'])
          end % end loop over time
         end % end loops over channels
       else   % rotate into RTZ !!!!
             
        for it=1:2:FT.nwin
         crap=zeros(3,1);X=zeros(Nx,1);Y=X;
         for ir=1:Rec(1).nr
         for ii=1:3;
               if mod(ir,Nx)==0;
                i1=floor(ir/Nx);i2=Nx;
                else
                 i1=floor(ir/Nx)+1;i2=mod(ir,Nx);
                end
          if ii <=2 
           crap(ii) = G_L(is,ir).U(ii,it,imode)+G_R(is,ir).U(ii,it,imode);
          else
           crap(ii) = G_R(is,ir).U(ii,it,imode);
          end
         end
         phi = Rec(ir).az / 180*pi;
         GG(i1,i2,1,it) = sin(phi)*crap(2)+cos(phi)*crap(1);
         GG(i1,i2,2,it) = cos(phi)*crap(2)-sin(phi)*crap(1);
         GG(i1,i2,3,it)= crap(3);
          X(i1)=Rec(ir).rg;
          Y(i2)=phi;
         end
        end
        [YY,XX]=meshgrid(Y,X);
        [AA,BB]=pol2cart(YY,XX);
        crap=zeros(Nx);
        for jj=1:3
          for it=1:2:1500%FT.nwin
           crap(1:Nx,1:Nx)=GG(1:Nx,1:Nx,jj,it);
              if it*FT.dt <=5
                max1=max(max(abs(crap)));
              end
           pcolor(AA,BB,crap)
           colormap(cc);shading interp;
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,0,'rp','Markersize',10,'Linewidth',3);axis square
            ylabel('Distance (km)');title(['channel ' char(chan(jj)) ' ' num2str(FT.tt(it)) ' seconds']);
            if it==1;pause(0.25);end
            [char(fdir) '/' char(chan(jj)) '_' int2str(it) '.jpg']
            caxis(max1/5*[-1 1])
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
      fdir = [char(fname) '/SOURCE_' int2str(Src(1).ns)];  
        mkdir(fdir)
        % arrange such that N = R, E=T and Z=Z
            for ii=1:3
            for it=1:2:FT.nwin
            time = FT.tt(it);
             for ir=1:Rec(1).nr
                 if mod(ir,Nx)==0;
                i1=floor(ir/Nx);i2=Nx;
             else
                 i1=floor(ir/Nx)+1;i2=mod(ir,Nx);
             end
             if ii <=2
              GG(i1,i2) = G_L(is,ir).U(ii,it,imode)+G_R(is,ir).U(ii,it,imode);
              GG(i1,i2) = GG(i1,i2)*sqrt(Rec(ir).rg);
                else
              GG(i1,i2) = G_R(is,ir).U(ii,it,imode);
              GG(i1,i2) = GG(i1,i2)*sqrt(Rec(ir).rg);
                end
              X(i1)=Rec(ir).rg;
              Y(i2)=Rec(ir).H;
            end
            [YY,XX]=meshgrid(Y,X);
            h1=pcolor(XX,2*YY,GG);axis ij;hold on
            if it*FT.dt <=5
                max1=max(max(abs(GG)));
            end
            colormap(cc);shading interp;axis equal;axis image
            %[h h]=contour(XX,YY,GG,[-1 0 1]/(10*max1),'k');
            set(gca,'Fontsize',14);xlabel('Distance (km)');hold on
            plot(0,Src(is).H,'rp','Markersize',10,'Linewidth',3);
            ylabel('Depth (km)');title(['channel ' char(chan(ii)) ' ' num2str(FT.tt(it)) ' seconds'])
            caxis(max1/10*[-1 1]);xlim([min(X) max(X)])
            [hh]=get(gca,'YTick');
            set(gca,'YTickLabel',num2str(hh))
            if it==1;pause(0.25);end
            [char(fdir) '/R_Z_' char(chan(ii)) '_' int2str(it) '.jpg']
            print('-djpeg','-r150',[char(fdir) '/' char(chan(ii)) '_' int2str(it) '.jpg'])
%             set(h,'LineStyle','none');
          %  delete(h);
            end % end loop over time
            end;end % end loops over channels
      
            
        end
    end
    
    
    
end
           