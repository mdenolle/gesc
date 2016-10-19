% Subroutine to solve eigenpropblem given a medium and frequency
% Marine Denolle (10/07/12)
 
 
function    [C,nerr] = eigSW(Cint,Dmax,omega,max_mode,poiss)


%% OUTPUTS
% C: structure that contains the medium properties per layer, the
% eigenfunctions, the dispersion curves etc etc
% nerr: if nerr = 0 ,  all good, if nerr = 1: error in solving eig pb
 
%% INPUTS
% Cint: structure that only contain the medium info describing from bottom
% to top
% Dmax: "halfspace" layer: as big as you want (pick at least 10xlongest
% wavelength), 
% omega: pulsation to solve eig for
% Nr : resolution number (number of points per wavelength)
% max_mode: maximum number of modes to look at
% if poiss == 1, lambda = mu, 0, lambda = rho.alpha^2 - 2*mu*beta^2

coo={'b.','r+','b.','r+','b.','r+','b.'};
C(1).Nl = Cint(1).Nl; 
C(1).omega=omega;
nerr = 0;
C(1).Ntot = 0;
if omega<1
    omega1 = 1;
else
    omega1 = omega;
end
 

for i=1:C(1).Nl
    
     % find appropriate number of points for each layer   
    if i==1 %halfspace
        C(i).H = Dmax-Cint(1).inter(1); % thickness
        C(i).N = 50; % fixed number of points
    else
        C(i).H = -Cint(1).inter(i)+Cint(1).inter(i-1);  % thickness      
        
        if i ~= C(1).Nl
            C(i).N = min([30 10+floor((1.1^log(omega1))*C(i).H)]);
        else
            C(i).N = min([40 floor((15*(1.1^log(omega1)))*(abs(log10(omega1+.1))+1.5))]);
            
        end
        
    end

    % now you have H and N, define Cheb D
    [C(i).D,x] = cheb(C(i).N-1);                            % differentiation matrix
    C(i).D=2/(C(i).H)*C(i).D;C(i).z=C(i).H/2*(x+1)';        % all that normalized - deep layer

    
    % fill in the zz, alpha, beta, rho and derived elastic parameters for
    % each layer and at each depth:
    for k=1:C(i).N
        if i==1  % for half space "layer"
           beta=Cint(1).beta(end);
           alpha=Cint(1).alpha(end);
           rho=Cint(1).rho(end);
        else
           b=Cint(i).betal;
           a=Cint(i).alphal;
           r=Cint(i).rhol;
           zz = Cint(1).inter(i)+C(i).z(k);  
           beta = interp1(Cint(i).zz,b,zz,'linear');
           rho = interp1(Cint(i).zz,r,zz,'linear');
           alpha = interp1(Cint(i).zz,a,zz,'linear');
        end     
        C(i).alpha(k)    = alpha;
        C(i).beta(k)     = beta;
        C(i).rho(k)      = rho;
        C(i).mu(k)       = rho*beta^2;
        C(i).lambda(k) =  rho.*alpha.^2 - 2*C(i).mu(k);
    end
    if poiss ==1 ; C(i).lambda = C(i).mu;C(i).beta=1/sqrt(3)*C(i).alpha;end
    C(i).lambdamu=C(i).lambda+2*C(i).mu;
   
    if i==1;
        C(1).beta1D = C(1).beta;C(1).alpha1D = C(1).alpha ; C(1).rho1D = C(1).rho;
        C(1).zz = C(1).z + Cint(1).inter(1);
    else
        C(1).beta1D  = [C(1).beta1D  , C(i).beta  ];
        C(1).alpha1D = [C(1).alpha1D , C(i).alpha ];
        C(1).rho1D   = [C(1).rho1D   , C(i).rho   ];
        C(1).zz      = [C(1).zz      , C(i).z+Cint(1).inter(i)];
    end
C(1).Ntot = C(1).Ntot+C(i).N;
end
C(1).mu1D = C(1).rho1D.*C(1).beta1D.^2;
C(1).lambda1D = C(1).rho1D.*C(1).alpha1D.^2 - 2*C(1).mu1D;
if poiss == 1;C(1).lambda1D = C(1).mu1D;C(1).beta1D=1/sqrt(3)*C(1).alpha1D;end
C(1).lambdamu1D = C(1).lambda1D + 2*C(1).mu1D;


%% solve eigenproblem
    % Rayleigh waves:
    tic;
        [A] = RayInnerMat(C,omega);
        [A,B] = RayBC(C,A); 
        [C,nerr] = EigRW(C,A,B,max_mode);
        if nerr ==1 ; return ; end

    % Love waves:
        [A] = LoveInnerMat(C,omega);
        [A,B] = LoveBC(C,A); 
        [C,nerr] = EigLW(C,A,B,max_mode);
        if nerr ==1 ; return ; end
        C(1).Nmode = min([max_mode length(C(1).kl) length(C(1).kr)]);

        
    %% get group velocity and integrals:
        [C]=get_integrals_sw(C);
        t1=toc;
        
        
