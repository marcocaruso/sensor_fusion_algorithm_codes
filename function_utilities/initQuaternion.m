function [z,qacc,qmag,l] = initQuaternion(a,m)
% Valenti 2015 - Keeping a Good Attitude: A Quaternion-Based Orientation
% Filter for IMUs and MARGs (Sensors, 2015)

% Implementation: Marco Caruso (Politecnico di Torino) 
% Date: 21/11/2019

% --------------- INPUT ---------------
% acc            = 1x3 (m/s^2)
% mag            = 1x3 (a.u.) normalized units

% --------------- OUTPUT ---------------
% z              = 1x4 [qw qx qy qz], the scalar part is at the beginning of the quaternion
% qacc           = 1x4 [qw qx qy qz], the scalar part is at the beginning of the quaternion
% qmag           = 1x4 [qw qx qy qz], the scalar part is at the beginning of the quaternion
% l              = 1x3 measured field rotated in the global horizontal plane
n=norm(a);
ax=a(1)/n; ay=a(2)/n; az=a(3)/n;
if az>=0
    qacc=[sqrt((az+1)/2) -ay/sqrt(2*(az+1)) ax/sqrt(2*(az+1)) 0]';
else
    qacc=[-ay/sqrt(2*(1-az)) sqrt((1-az)/2) 0 ax/sqrt(2*(1-az))]';
end

hx=m(1); hy=m(2); hz=m(3);
l=quatrotmatr(qacc)'*[hx;hy;hz];
T=l(1)^2+l(2)^2;
if l(1)>=0
    qmag=[sqrt(T+l(1)*sqrt(T))/sqrt(2*T);...
        0;0;...
        l(2)/(sqrt(2)*sqrt(T+l(1)*sqrt(T)))];
else
    qmag=[l(2)/(sqrt(2)*sqrt(T-l(1)*sqrt(T)));...
        0;0;...
        sqrt(T-l(1)*sqrt(T))/sqrt(2*T)];
end
z=quatmultiply(qacc',qmag');

end