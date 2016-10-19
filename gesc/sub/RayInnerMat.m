function [A] = RayInnerMat(C,omega)
Nl = C(1).Nl;   % number of layers


A=zeros(C(1).Ntot);
n = 0;
for i=1:Nl
    I=(n+1:n+4*C(i).N);
    A(I,I) = [zeros(C(i).N) , -diag(C(i).lambda./C(i).lambdamu)*C(i).D , diag(1./C(i).lambdamu) , zeros(C(i).N) ; ...
                    C(i).D , zeros(C(i).N) , zeros(C(i).N) , -diag(1./C(i).mu); ... 
                    diag(C(i).rho)*(omega*omega') , zeros(C(i).N),zeros(C(i).N),C(i).D ; ....
                    zeros(C(i).N),-diag(C(i).rho)*(omega*omega') - C(i).D*diag( 4*C(i).mu.*(C(i).lambda+C(i).mu)./C(i).lambdamu)*C(i).D ,...
                       -C(i).D*diag(C(i).lambda./C(i).lambdamu),zeros(C(i).N)]; % bottom layer            
   n = n + 4*C(i).N;
end

