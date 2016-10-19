function [C,nerr]=EigRW(C,A,B,Nmode)
nerr = 0;
Nl = C(1).Nl;

%% solve GEP
% [Ll,Lr] = RaySel(C);

A(isnan(A)) = 0;
A(isinf(A)) = 0;
B(isnan(B)) = 0;
B(isinf(B)) = 0;
Anew=A;
Bnew=B;

% Anew = (Ll*A*Lr);

% Bnew = (Ll*B*Lr);
% if min(C(1).mu1D)==0 % if there is a water layer at the top
% %     % impose continuity condition for r2
% disp('got water')
%     n=C(end).N;
%    Anew(end-2*n+1,:)=0;
%    Anew(end-2*n+1,end-2*n-C(end-1).N)=-1;
%    Anew(end-2*n+1,end-2*n)=1;
%    Bnew(end-2*n+1,:)=0;
%    
%    % R2=0
%    Anew(end-C(end).N,:)=0;
% Bnew(end-C(end).N,:)=0;
% Anew(end-C(end).N,end-3*C(end).N+1:end-2*C(end).N) = -C(end).lambdamu(end).*C(end).D(end,:); % )
% Bnew(end-C(end).N,end-3*C(end).N                 ) = +C(end).lambda  (end);
% 
% 
% end

[X,k] = eig(Anew,Bnew);
%[X,k] = eig(A,B);

k = diag(k);
ik = find(real(k)./imag(k)>=1E15&real(k)~=Inf&real(k)>0&C(1).omega./real(k)>0.1) ;
k = k(ik) ; X = X(:,ik);
[k,ii] = sort(k,'descend'); X = X(:,ii);
% figure(301)
% plot(C(1).omega/2/pi./k,'o');grid on
% pause
disp(['# of Rayleigh modes found: ' int2str(length(k))])
if isempty(k) 
    nerr = 1;return
end
%disp('magic')
%size(X)
% X = Lr*X;
%size(X)
%figure; plot(X(:,1))
%disp('magic')

%figure; imagesc(Lr); 

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


%% eigenfunctions:
C(1).ux = zeros(C(1).Ntot,1);
C(1).uz = zeros(C(1).Ntot,1);
C(1).r3 = zeros(C(1).Ntot,1);
C(1).r4 = zeros(C(1).Ntot,1);
C(1).szz= zeros(C(1).Ntot,1);

for ik=1:min([length(k) Nmode])
    n=0;
    toto1=X(0*C(1).N+1:1*C(1).N,ik)';
    toto2=X(1*C(1).N+1:2*C(1).N,ik)';
    toto3=X(2*C(1).N+1:3*C(1).N,ik)';
    toto4=X(3*C(1).N+1:4*C(1).N,ik)';
    for i=2:Nl
      n = n + 4*C(i-1).N;
      toto1 = [toto1 , X(n+0*C(i).N+1:n+1*C(i).N,ik)'];
      toto2 = [toto2 , X(n+1*C(i).N+1:n+2*C(i).N,ik)'];
      toto3 = [toto3 , X(n+2*C(i).N+1:n+3*C(i).N,ik)'];
      toto4 = [toto4 , X(n+3*C(i).N+1:n+4*C(i).N,ik)'];
    end
    
    normind=C(1).Ntot;
%     for ii=Nl:-1:1; if (C(1).solid(ii)~=1); normind=normind-C(ii).N; break ; end; end;
    %[null,normind] = min(C.solid1D); normind = normind-1;
    %if (min(C(1).solid1D)==1); normind=length(C(1).solid1D); end;
    %if (max(C(1).solid1D)==0); normind=length(C(1).solid1D); end;
%     normalize=max(abs(toto2))*sign(toto2(end));%(normind);
%     if ik==1;normalize=toto2(end);end
%     if (normalize==0); normalize = max(toto2); end;
    normalize=toto2(end);
    C(1).ux(:,ik) = toto1/normalize;
    C(1).uz(:,ik) = toto2/normalize;
    C(1).r3(:,ik) = toto3/normalize;
    C(1).r4(:,ik) = toto4/normalize;
%     if ik==1&& min(toto2)*max(toto2)<0&&min(toto2)<-0.1;
%             disp('modes flipped')
%     end

%     CC.ux(:,ik) = toto1;
%     CC.uz(:,ik) = toto2;
%     CC.r3(:,ik) = toto3./CC.ux(end,ik);
%     CC.r4(:,ik) = toto4./CC.ux(end,ik);
%     CC.uz(:,ik) = CC.uz(:,ik)./CC.ux(end,ik);CC.ux(:,ik)=CC.ux(:,ik)./CC.ux(end,ik); 

    C(1).kr(ik)=k(ik);

    C(1).szz(:,ik) = 4*(C(1).mu1D.*(C(1).lambda1D+C(1).mu1D)./C(1).lambdamu1D).*(C(1).BigD*C(1).uz(:,ik))' ...
        + C(1).lambda1D./C(1).lambdamu1D.*C(1).r3(:,ik)';
    
    
end
