% Subroutine to solve eigenpropblem given a medium and frequency
% Marine Denolle (04/10/14)
 
 
function    bigC = gesc(MED_GESC,FT,NR,max_mode,freq)

%% OUTPUTS
% bigC: structure that contains the medium properties per layer, the
% eigenfunctions, the dispersion curves etc etc
 
%% INPUTS
% MED:medium structure
% FT: time/frequency structure
% NR: structure of variation of number of collocation points with
% wavelength (contains for upper layer
% NR : number of points scheme: see Section 2.3 on resolution 
%     NR.typ = 'cst','lin' or 'exp' for type of variation of #points with
%     thickness
%     NR.rate = 'cst' (constant slope), 'rec' (recursive and depends on
%     previous variations of the eigen functions)
%     NR.lambda = depths at which change scheme to lower resolution number 
%     NR.amp(1:length(NR.lambda)) = slope value (see Figure XX from manual)
%     NR.N1(1:length(NR.lambda)) = lower bound # points
%     NR.N2(1:length(NR.lambda)) = upper bound # points

% max_mode: maximum number of modes to look at
% freq(optional): user may ask to output eigenfunctions at frequencies stored in
% the freq vector.
% bigC(1:min(length(FT.omega),floor(FT.Fmax/FT.df+1)))=struct('a',[]);
bigC(1:min(length(FT.omega),floor(FT.Fmax/min(FT.df)+1)))=struct('a',[]);
% calculate for each frequency, up to Nyquist or Fmax.
if length(FT.df)==1&&floor(FT.Fmax/FT.df)>length(FT.omega);
    disp('!! Choose Fmax < Fnyq !!')
    return
end
nerr=1;
 
    
for ifreq = 1:min(length(FT.omega),floor(FT.Fmax/min(FT.df)+1))
    clear C
    if FT.omega(ifreq)==0;continue;end
tic;
MED=MED_GESC;

% If most of the eigenfunction variations is trapped within the upper
% layer, we split that layer to reduce the number of points required to
% solve the strong variations in the layer. Dichotomy process!

if ifreq==1 || nerr==1;
    C1(1).cl(1)=max(MED(1).beta);
    C1(1).cr(1)=max(MED(1).beta);
    C1(1).Nl=MED(1).Nl;
    HH =  max(MED_GESC(1).inter);
else
ikk=find(C1(1).uz(:,1)<=1E-2*abs(C1(1).uz(end,1)));
HH=C1(1).zz(ikk(end));
end
if ifreq>1 && HH<MED_GESC(1).inter(end-1)

   
        
if strcmp(MED_GESC(1).typ,'cst')  % if this a layered medium with homoheneous layers
    for i=1:length(MED_GESC(1).inter)
        [~,ikk]=min(abs(MED_GESC(1).z-MED_GESC(1).inter(i)));
        aa(i) = MED_GESC(1).alpha(ikk);
        bb(i) = MED_GESC(1).beta(ikk);
        rr(i) = MED_GESC(1).rho(ikk);
    end
    aa(i+1) = aa(i);bb(i+1) = bb(i);rr(i+1) = rr(i);
    aa=fliplr(aa);bb=fliplr(bb);rr=fliplr(rr);
    MED(1).inter = [0 HH fliplr(MED_GESC(1).inter(1:end-1))];
    MED= make_layers(MED(1).inter,MED_GESC(1).Dmax,aa,bb,rr,MED_GESC(1).typ);

elseif strcmp(MED_GESC(1).typ,'lin')  % if you have gradients
    i=MED_GESC(1).Nl;
    MED(1).inter=[MED_GESC(1).inter(1:end-1) HH 0];
    MED(1).Nl=length(MED(1).inter);
    for i=MED(1).Nl-1:MED(1).Nl
        MED(i).zz=linspace(MED(1).inter(i-1),MED(1).inter(i),30);
        MED(i).alphal = interp1(MED_GESC(1).z,MED_GESC(1).alpha,MED(i).zz,'linear');
        MED(i).betal  = interp1(MED_GESC(1).z,MED_GESC(1).beta, MED(i).zz,'linear');
        MED(i).rhol   = interp1(MED_GESC(1).z,MED_GESC(1).rho,  MED(i).zz,'linear');

    end
end
      disp(['Added a layer because first wavelength is ' num2str(HH)  ...
        'km compared to ' num2str(MED(1).inter(end-1)) 'km'])
end
C(1).Nl = MED(1).Nl;
C(1).omega=FT.omega(ifreq);
C(1).Ntot = 0;
C(1).solid=MED(1).solid;

T=2*pi/C(1).omega;
disp(['Period: ' num2str(T) ' s and Frequency ' num2str(C(1).omega/(2*pi)) ' Hz'])

    
for i=1:C(1).Nl
    
     % find appropriate number of points for each layer   
    if i==1 %halfspace
        C(i).H = MED(1).Dmax-MED(1).inter(1); % thickness
    else
        C(i).H = -MED(1).inter(i)+MED(1).inter(i-1);  % thickness
    end

    
        if strcmp(NR(1).typ,'typ')==1
        LL = C(i).H/(C1(1).cl(1)*T); % wavelength
      [~,ii]=find(MED(1).inter(i)/(C1(1).cl(1)*T)<=[NR.lambda 1E4]);
      [~,kk]=min(abs(LL-NR(1).x));
      C(i).N = floor(NR(ii(1)).npts(kk));
        end
        if strcmp(NR.typ,'rec')==1
            mgrd=1;
            if ifreq>1
            % find peak gradient in that layer
            if (C(1).Nl==C1(1).Nl || i < C1(1).Nl) && nerr==0
             for imode=1:1
                n1=(C1(i).nn(1)+1:C1(i).nn(2));
                u1=C1(1).ux(n1,imode);u2=C1(1).uy(n1,imode);u3=C1(1).uz(n1,imode);
                dr1 = max(abs(C1(i).D*u1));
                dl1 = max(abs(C1(i).D*u2));
                dr2 = max(abs(C1(i).D*u3));
               grd(imode) = mean([dr1 dr2 dl1]);  
            end
            mgrd=mean(grd);
            end
            end
            C(i).N =min(max(floor(NR.N1*mgrd*C(i).H),NR.N1),NR.N2);
        end
      
    % now you have H and N, define Cheb D
    [C(i).D,x] = cheb(C(i).N-1);                            % differentiation matrix
    C(i).D=2/(C(i).H)*C(i).D;C(i).z=C(i).H/2*(x+1)';        % all that normalized - deep layer

    
    % fill in the zz, alpha, beta, rho and derived elastic parameters for
    % each layer and at each depth:
    for k=1:C(i).N
          [~,ii]=unique(MED(i).zz,'stable');
        if i==1  % for half space "layer"
           beta=max(MED(1).beta);
           alpha=max(MED(1).alpha);
           rho=max(MED(1).rho);
        else
           b=MED(i).betal(ii);
           a=MED(i).alphal(ii);
           r=MED(i).rhol(ii);
           zz = MED(1).inter(i)+C(i).z(k);  
           beta = interp1(MED(i).zz(ii),b,zz,'linear');
           rho = interp1(MED(i).zz(ii),r,zz,'linear');
           alpha = interp1(MED(i).zz(ii),a,zz,'linear');
        end
        C(i).alpha(k)    = alpha;
        C(i).beta(k)     = beta;
        C(i).rho(k)      = rho;
        C(i).mu(k)       = rho*beta^2;
        C(i).lambda(k)   = rho.*alpha.^2 - 2*C(i).mu(k);
            
    end
    C(i).lambdamu=C(i).lambda+2*C(i).mu;
    if i==1;
        C(1).beta1D = C(1).beta;C(1).alpha1D = C(1).alpha ; C(1).rho1D = C(1).rho;
        C(1).zz = C(1).z + MED(1).inter(1);
    else
        C(1).beta1D  = [C(1).beta1D  , C(i).beta  ];
        C(1).alpha1D = [C(1).alpha1D , C(i).alpha ];
        C(1).rho1D   = [C(1).rho1D   , C(i).rho   ];
        C(1).zz      = [C(1).zz      , C(i).z+MED(1).inter(i)];
    end
C(1).Ntot = C(1).Ntot+C(i).N;
Ntot(ifreq)=C(1).Ntot;
end

C(1).mu1D = C(1).rho1D.*C(1).beta1D.^2;
C(1).lambda1D = C(1).rho1D.*C(1).alpha1D.^2 - 2*C(1).mu1D;
% if poiss == 1;C(1).lambda1D = C(1).mu1D;C(1).beta1D=1/sqrt(3)*C(1).alpha1D;end
C(1).lambdamu1D = C(1).lambda1D + 2*C(1).mu1D;

if ifreq>1 && C(1).Ntot > 500;return;end
%% solve eigenproblem
% Rayleigh waves:
    [A] = RayInnerMat(C,C(1).omega);
    [A,B] = RayBC(C,A);
    [C,nerr] = EigRW(C,A,B,max_mode);
    if (nerr==1);continue;end
% Love waves:
    [A] = LoveInnerMat(C,C(1).omega);
    [A,B] = LoveBC(C,A); 
    [C,nerr] = EigLW(C,A,B,max_mode);
    if (nerr==1);continue;end
    C(1).Nmode = min([max_mode length(C(1).kl) length(C(1).kr)]);


%% get group velocity and integrals:
    [C]=get_integrals_sw(C);
    t1=toc;
    tt(ifreq)=t1;
       
    C1=C;
    bigC(ifreq) = struct('a',C);
    if nargin>4
       [aa,ikk]=min(abs(freq-C(1).omega/(2*pi)));
       if aa>max(FT.df)/2;continue;end
       figure(3)
       for imode=1:C(1).Nmode
        plot(C(1).ux(:,imode),(C(1).zz),'b+-',C(1).uy(:,imode),(C(1).zz),'k--',...
            C(1).uz(:,imode),(C(1).zz),'r-.','Markersize',8,'Linewidth',2);hold on;grid on
        legend('r_1','l_1','r_2')
       end
       axis ij
       
       title(['Eigenfunctions at T = ' num2str(2*pi/C(1).omega) ' s' ])
       set(gca,'Fontsize',14);ylabel('Deph (km)');axis ij
       ylim([0 50])
       print('-dpsc',['eigenfunctions_T' num2str(floor(2*pi/C(1).omega)) '.ps']);
%        close all
    end
    [~,ib]=unique(C(1).zz);
    uyi(ifreq) = interp1(C(1).zz(ib),C(1).uy(ib),5,'linear');
%     disp(uyi(ifreq))
%     figure(100)
%     subplot(211)
%     plot(C(1).zz(ib),C(1).uy(ib));xlim([0 10]);grid on
%     subplot(212)
%     semilogx(FT.omega(1:ifreq)/2/pi,uyi(1:ifreq),'r+');grid on;pause(0.01)
end

