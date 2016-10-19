%% make collocation sampling rule: Number of points / wavelength.
function Nr = make_npt_colloc_lin(lambda,N1,N2,a1,a2)

%% inputs:
% lambda: vector of depth segment to attribute different scaling
% N1, N2: lower and upper bounds of number of points (size lambda-1)
% a1, a2: start scaling relation between a1 x lamba and a2 x lambda
%% output:
% Nr: number of points
Nr=struct;
Nr(1).lambda=lambda;
Nr(1).a1=a1;Nr(1).a2=a2;
Nr(1).a1=N1;Nr(1).a2=N2;
if length(N1)~=length(lambda);disp('check resolution input');return;end
x=logspace(-3,3,100); % vector of thickness / wavelength

Nr(1).typ='lin';
for i=1:length(lambda)
    Nr(i).lambda=lambda(i);
    Nr(i).x=x;
    Nr(i).npts(1:length(x)) =  min(max((N2(i)-N1(i))/(a2(i)-a1(i))*x,N1(i)),N2(i));
end
Nr(i+1).lambda=1E6;
Nr(i+1).npts = N1(i)*ones(length(x),1); 



return