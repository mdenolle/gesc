function [A,B] = LoveBC(C,A)



Nl = length(C);   % number of layers
B=eye(size(A));   % 

% Bottom Boundary Condition (no displacement => remove an equation of motion)
A(1,:)=0;
A(1,1) = 1 ; B(1,1) = 0;
% Top Boundary Condition (no-traction => remove another equation of motion)
A(end-C(Nl).N,:)=0;
A(end-C(Nl).N,end-2*C(Nl).N+1:end-C(Nl).N) = C(Nl).D(C(Nl).N,:) ; B(end-C(Nl).N,:)=0;

n = 0;
for i=1:Nl-1      % foreach interface
    n = n + 2*C(i).N;
    % stress
    A(n-C(i).N,:)=0;B(n-C(i).N,:)=0;
    A(n-C(i).N,n-2*C(i).N+1:n-C(i).N)=-C(i).mu(C(i).N)*C(i).D(C(i).N,:) ;
    A(n-C(i).N,n+1:n+C(i+1).N) = C(i+1).mu(1)*C(i+1).D(1,:) ;
    % displacement
    A(n+1,:)=0;
    A(n+1,n-C(i).N)=-1;A(n+1,n+1)=1;B(n+1,:)=0;
   
end



% 
%     A(1,:) = 0; A(1,1) = 1;B(1,1)=0;
%     % zero traction at surface: Duy/dz = 0
%     A(2*N1+N2,:) = 0; A(2*N1+N2,2*N1+1:2*N1+N2) = C(2).D(N2,:);B(2*N1+N2,2*N1+N2)=0;
%     % continuity displacement at interface:
%     A(2*N1+1,:)=0;A(2*N1+1,N1)=-1;A(2*N1+1,2*N1+1)=1;B(2*N1+1,:)=0;
%     % continuity of stress
%     A(N1,:)=0;A(N1,1:N1) = -C(1).mu(N1)*C(1).D(N1,:);A(N1,2*N1+1:2*N1+N2) = C(2).mu(1)*C(2).D(1,:);B(N1,N1)=0;
%     