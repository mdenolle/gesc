% Subroutine to get the surface waves integrals products and group velocity
% from  density and elastic profiles and Eigenfunctions obtained from
% Chebyshev collocation technique 
% Marine Denolle (01/27/11)



function [C]=get_integrals_sw(C)


bigD = C(1).BigD;
W  =C(1).W;

rho=C(1).rho1D;
mu=C(1).mu1D;
lambda=C(1).lambda1D;
lambdamu=C(1).lambdamu1D;

    
ux2=zeros(length(C(1).zz),C(1).Nmode);
uy2=ux2;uz2=ux2;C(1).Dux=ux2;C(1).Duy=ux2;C(1).Duz=ux2;
C(1).sigmazz=ux2;C(1).sigmazz2=ux2;C(1).r32=ux2;
for i=1:C(1).Nmode
        D1ux = bigD*C(1).ux(:,i);
        D1uy = bigD*C(1).uy(:,i);
        D1uz = bigD*C(1).uz(:,i);
        
        ux2(:,i) = C(1).ux(:,i).^2;
        uy2(:,i) = C(1).uy(:,i).^2;
        uz2(:,i) = C(1).uz(:,i).^2;

       % integrals:
       %Love:
        C(1).Il(1,i) = 1/2*dot(W,rho.*uy2(:,i)');
        C(1).Il(2,i) = 1/2*dot(W,mu.*uy2(:,i)');
        C(1).Il(3,i) = 1/2*dot(W,mu.*(D1uy.^2)');
        C(1).Ul(i) = C(1).Il(2,i)*C(1).kl(i)/(C(1).Il(1,i)*C(1).omega);
        C(1).cl(i) = C(1).omega/C(1).kl(i);
        
        % Rayleigh:
        C(1).Ir(1,i) = 1/2*dot(W,rho.*(ux2(:,i)'+uz2(:,i)'));
        C(1).Ir(2,i) = 1/2*dot(W,lambdamu.*ux2(:,i)'+mu.*uz2(:,i)');
        C(1).Ir(3,i) = dot(W,lambda.*(C(1).ux(:,i).*D1uz)' - mu.*(C(1).uz(:,i).*D1ux)');
        C(1).Ir(4,i) = 1/2*dot(W,lambdamu.*(D1uz.^2)'+mu.*(D1ux.^2)');
        
        C(1).cr(i) = C(1).omega/C(1).kr(i);c=C(1).cr(i);
        C(1).Ur(i) = (C(1).Ir(2,i)+C(1).Ir(3,i)/(2*C(1).kr(i)))/(c*C(1).Ir(1,i));

        C(1).Dux(:,i) = D1ux;       
        C(1).Duy(:,i) = D1uy;
        C(1).Duz(:,i) = D1uz;
end

%         disp('checking sizes')
%         disp([size(C(1).Dux) size(C(1).ux)])
