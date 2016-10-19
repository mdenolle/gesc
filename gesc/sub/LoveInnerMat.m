function [A] = LoveInnerMat(C,omega)


Nl = C(1).Nl;   % number of layers

ntot=0;
for i=1:Nl
    ntot = ntot + 2*C(i).N;
end
A=zeros(ntot);


n = 0;
for i=1:Nl
   I = (n+1:n+2*C(i).N);
   A(I,I) =  [zeros(C(i).N) , diag(1./C(i).mu) ; ...
                    C(i).D*diag(C(i).mu)*C(i).D + diag(C(i).rho)*omega^2 , zeros(C(i).N)];
 
    % fill matrix
   %  A(n+2:n+C(i).N-1,n+1+C(i).N:n+2*C(i).N) =  [zeros(C(i).N-2,1) , diag(1./C(i).mu(2:C(i).N-1)) , zeros(C(i).N-2,1) ];
%     A(n+2:n+C(i).N-1,n+1+C(i).N:n+2*C(i).N) =  diag(1./C(i).mu(2:C(i).N-1));
   % A(n+C(i).N+1:n+2*C(i).N,n+1:n+C(i).N) = omega^2*diag(C(i).rho(1:C(i).N)) + C(i).D*diag(C(i).mu(1:C(i).N))*C(i).D;
   % A(n+C(i).N+1,:)=0;A(n+2*C(i).N,:)=0;
    
   n = n + 2*C(i).N;
end
