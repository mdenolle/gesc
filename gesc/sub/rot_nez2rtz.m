 % function to rotate green tensor from NEZ to RTZ

function G2 = rot_nez2rtz(G,az)

% NEZ => z downward direct system
% RTZ => Z downward direct system


phi = az*pi/180;
a = size(G);

R = [cos(phi) sin(phi) 0 ; ...
    -sin(phi) cos(phi) 0 ;...
    0 0 1 ] ;



if a(1)>1 && a(2)>1
    if min(a) == 3 && max(a) > 3 
    G2(3,:)=G(3,:);
    G2(1,:) = sin(phi)*G(2,:)+cos(phi)*G(1,:);
    G2(2,:) = cos(phi)*G(2,:)-sin(phi)*G(1,:);
    else
    G2 = R*G*R';
    end
end
  
  % Z = > 3
  % R = > 1
  % T = > 2
  
  %[zz zr zt ; rz rr rt; tz tr tt]
