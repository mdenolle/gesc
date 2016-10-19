function [A,B] = RayBC(C,A)
Nl = C(1).Nl;   % number of layers
B  = eye(size(A));

% Bottom BC's (remove two above)

% Uz =0 bottom
A(1+1*C(1).N,:) = 0;
A(1+1*C(1).N,1) = 1;
B(1+1*C(1).N,:) = 0;

% Ux = 0
A(1+0*C(1).N,         :) = 0;
A(1+0*C(1).N,1+1*C(1).N) = 1;
B(1+0*C(1).N,         :) = 0;



% for each interface (replace one above and one below the interface)
%disp('series 1')
n = 0;
for ii=1:Nl-1;
    n = n + 4*C(ii).N;
    % Continuity conditions
    % stress sigma_zz
    A(n-C(ii).N,:) = 0;
    B(n-C(ii).N,:) = 0;
    A(n-C(ii).N,n-3*C(ii  ).N+1:n-2*C(ii  ).N) = -C(ii  ).lambdamu(end)*C(ii  ).D(end,:); 
    B(n-C(ii).N,n-3*C(ii  ).N                ) = +C(ii  ).lambda  (end);
    A(n-C(ii).N,n+  C(ii+1).N+1:n+2*C(ii+1).N) = +C(ii+1).lambdamu(1  )*C(ii+1).D(1  ,:);
    B(n-C(ii).N,n+           1               ) = -C(ii+1).lambda  (1  );
    
    % Continuity conditions on sigma_xz (R4)
    A(n,:) =  0;
    A(n,n) = -1;
    A(n,n+3*C(ii+1).N+1) = 1;
    B(n,:) =  0;    % 
 
    
    % continuous r1 (ux)
    A(n+1,            :) =  0;
    A(n+1,n-3*C(ii  ).N) = -1;
    A(n+1,n+1) =  1;
    B(n+1          ,  :) =  0; 
    
    % Continuity conditions on uz (r2)
    A(n+1+C(ii+1).N,          :) =  0;
    A(n+1+C(ii+1).N,n-2*C(ii).N) = -1;
    A(n+1+C(ii+1).N,n+1+C(ii+1).N) =  1;
    B(n+1+C(ii+1).N,          :) =  0; 

        
    % if the interface is water-solid, free shear condition
    if max(C(ii+1).mu)==0
        % Continuity conditions on sigma_xz (R4)
%         A(n,n+3*C(ii+1).N+1) = 0;
    end
    
end

% BC's and continuity condition on sigma_xz (R4) and ux (r2)
% %disp('series 2')
% n = 0;
% for ii=1:Nl-1;
%     n = n + 4*C(ii).N;
%     
%     % Interface between two solids (replace one above and one below the interface)
%     if(C(1).solid(ii) == C(1).solid(ii+1) && C(1).solid(ii) == 1);
% 
%         % Continuity conditions on sigma_xz (R4)
%         A(n,:) =  0;
%         A(n,n) = -1;
%         A(n,n+3*C(ii+1).N+1) = 1;
%         B(n,:) =  0;    % 
%         
%         % Continuity conditions on uz (r2)
%         A(n+1+C(ii+1).N,          :) =  0;
%         A(n+1+C(ii+1).N,n-2*C(ii).N) = -1;
%         A(n+1+C(ii+1).N,n+1        ) =  1;
%         B(n+1+C(ii+1).N,          :) =  0; 
%         
%         
% %         % Continuity conditions on ux (r1)
%         A(n+1       ,          :) =  0;
%         A(n+1       ,n-3*C(ii).N) = -1;
%         A(n+1       ,n+1        ) =  1;
%         B(n+1       ,          :) =  0; 
%         continue;
%     end;
%     
%     % sigma_xz (R4=0)
%     A(n,:) = 0;
%     A(n,n) = 1;
%     B(n,:) = 0;         % sigma_xz==0   
%     break;
% end

% Top BC's 
% free traction: impose on r3 the condition sigma_zz=0
A(end-C(end).N,:)=0;
B(end-C(end).N,:)=0;
A(end-C(end).N,end-3*C(end).N+1:end-2*C(end).N) = -C(end).lambdamu(end).*C(end).D(end,:); % )
B(end-C(end).N,end-3*C(end).N                 ) = +C(end).lambda  (end);
 % sigma_xz (R4=0) sigma_xz==0
%  disp(C(end).solid)
%  if (C(end).solid==1)
A(end,:) = 0;
A(end,end) = 1;
B(end,:) = 0;
%  end


