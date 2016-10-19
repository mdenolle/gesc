function [C,nerr]=EigLW(C,A,B,Nmode)
nerr = 0;
Nl=C(1).Nl;

% solve GEP
% [L] = LoveSel(C);

%size(L)
%size(A)

A(isnan(A)) = 0;
A(isinf(A)) = 0;
B(isnan(B)) = 0;
B(isinf(B)) = 0;

Anew=A;
Bnew=B;
% Anew = L*A*L';
% Bnew = L*B*L';

[X,k] = eig(Anew,Bnew);
%[X,k] = eig(A,B);

k = diag(k);
ik = find(imag(k)==0&real(k)>0&real(k)~=Inf) ; k = k(ik) ; X = X(:,ik);
[k,ii] = sort(k,'descend');X = X(:,ii);
% figure(100)
% plot(C(1).omega./k)
disp(['# of Love modes found: ' int2str(length(k))])
if isempty(k) 
    nerr = 1;return
end
% X = L'*X;

n1=0;W=[];bigD=0;
for i=1:C(1).Nl
    n2 = n1 + C(i).N;
    bigD(n1+1:n2,n1+1:n2)=C(i).D;
    [~,W(n1+1:n2)]=clencurt(C(i).N-1);
    W(n1+1:n2) = C(i).H/2*W(n1+1:n2);
    C(i).nn=[n1 n2];
    n1=n2;
end

C(1).BigD = bigD;
C(1).W = W;

bigD = C(1).BigD;
W  =C(1).W;

% eigenfunctions
C(1).uy  = zeros(C(1).Ntot,1);
C(1).syz = zeros(C(1).Ntot,1);

toto=[];
for ik=1:min([length(k) Nmode])
    n    = 0;
    toto = X(1:C(1).N,ik)';
    for i=2:Nl
      n    = n + 2*C(i-1).N;
      toto = [toto , X(n+1:n+C(i).N,ik)'];
    end
    
    normind=C(1).Ntot;
%     for ii=Nl:-1:1; if (C(1).solid(ii)~=1); normind=normind-C(ii).N; break ; end; end;
%     normalize=max(abs(toto))*sign(toto(end));%(normind);
%     if ik==1;normalize=toto(end);end
%     if (normalize==0); normalize = max(toto); end;
    normalize=toto(end);
    
    C(1).uy(:,ik) = toto./normalize;
    C(1).kl(ik)=k(ik);
    C(1).syz(:,ik) = (C(1).mu1D.*(C(1).BigD*C(1).uy(:,ik))')';     
end
