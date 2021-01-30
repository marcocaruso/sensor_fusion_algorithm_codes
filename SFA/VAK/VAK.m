function [qvak] = VAK(XX, t, stdGyro, stdAcc, stdMag)

% "A Linear Kalman Filter for MARG Orientation Estimation Using the Algebraic Quaternion Algorithm" by (Roberto G. Valenti, Ivan Dryanovski and Jizhong Xiao)


% --------------- INPUT ---------------
% XX             = (:,2:4) Accelerometer (m/s^2), (:,5:7) Gyroscope (rad/s), (:,8:10) Magnetometer (a.u.) normalized units
% t              = Nx1 time vector (s)
% stdAcc         = 1x1 (m/s^2) inverse weight to assign to the accelerometer measurements
% stdGyro        = 1x1 (rad/s) inverse weight to assign to the gyroscope measurements
% stdMag         = 1x1 (a.u.)  inverse weight to assign to the magnetometer measurements

% --------------- OUTPUT ---------------
% qvak           = Nx4 [qw qx qy qz]

% Implementation by Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 19/05/2020



%Accelerometer
ax = -XX(:,2); ay = -XX(:,3); az = -XX(:,4);

%Gyroscope
wx = XX(:,5); wy = XX(:,6); wz = XX(:,7);

%Magnetometer
hx = XX(:,8); hy = XX(:,9); hz = XX(:,10);

sigmaU = [stdAcc^2*eye(3) zeros(3)
    zeros(3) stdMag^2*eye(3)];

dt = mean(t);

x_priori = zeros(4,length(t));
x_post = zeros(4,length(t));
x_post(:,1) = MeasQuaternion([ax(1) ay(1) az(1)],[hx(1) hy(1) hz(1)]);

CSI = [x_post(2,1) x_post(3,1) x_post(4,1)
    -x_post(1,1) -x_post(4,1) -x_post(3,1)
    x_post(3,1) -x_post(1,1) -x_post(2,1)
    -x_post(3,1) x_post(2,1) -x_post(1,1)];

Q = (dt/2)^2*CSI*stdGyro*CSI';
P_post = Q;

for i=1:length(t)-1
    dt = t(i+1) - t(i);
    
    omega = 0.5*[0 wx(i) wy(i) wz(i)
        -wx(i) 0 wz(i) -wy(i)
        -wy(i) -wz(i) 0 wx(i)
        -wz(i) wy(i) -wx(i) 0];
    
    phi = eye(4) + omega*dt;
    
    x_priori(:,i) = phi*x_post(:,i);
    
    CSI = [x_post(2,i) x_post(3,i) x_post(4,i)
        -x_post(1,i) -x_post(4,i) -x_post(3,i)
        x_post(3,i) -x_post(1,i) -x_post(2,i)
        -x_post(3,i) x_post(2,i) -x_post(1,i)];
    
    
    Q = (dt/2)^2*CSI*stdGyro*CSI';
    P_priori = phi*P_post*phi'+Q;
    
    [z,qacc,qmag,l] = MeasQuaternion([ax(i) ay(i) az(i)],...
        [hx(i) hy(i) hz(i)]);
    if isrow(z)
        z = z';
    end

    J = CovEstim(qacc,qmag,[ax(i) ay(i) az(i)],l,[hx(i) hy(i) hz(i)]);
    R = J*sigmaU*J';

    %Compute the Kalman Gain
    warning off
    H = eye(4);
    K = P_priori*H'*(H*P_priori*H'+R)^-1;
    warning on
   
    if dot(x_priori(:,i),z)<0
        z = -z;
    end
    
    % Compute the a posteriori state estimate
    x_post(:,i+1) = x_priori(:,i)+K*(z-x_priori(:,i));
    
    % Normalize the quaternion
    x_post(1:4,i+1) = x_post(1:4,i+1)/norm(x_post(1:4,i+1));
    
    %Compute the a posteriori covariance matrix
    P_post = P_priori-K*H*P_priori;
end

qvak = x_post';
end

function [z,qacc,qmag,l] = MeasQuaternion(a,m)
n = norm(a);
ax = a(1)/n; ay = a(2)/n; az = a(3)/n;
if az>=0
    qacc = [sqrt((az+1)/2) -ay/sqrt(2*(az+1)) ax/sqrt(2*(az+1)) 0]';
else
    qacc = [-ay/sqrt(2*(1-az)) sqrt((1-az)/2) 0 ax/sqrt(2*(1-az))]';
end
clear n

n = norm(m);
hx = m(1)/n; hy = m(2)/n; hz = m(3)/n;
l = quatrotmatr(qacc)'*[hx;hy;hz];
T = l(1)^2+l(2)^2;
if l(1) >= 0
    qmag = [sqrt(T+l(1)*sqrt(T))/sqrt(2*T);...
        0;0;...
        l(2)/(sqrt(2)*sqrt(T+l(1)*sqrt(T)))];
else
    qmag = [l(2)/(sqrt(2)*sqrt(T-l(1)*sqrt(T)));...
        0;0;...
        sqrt(T-l(1)*sqrt(T))/sqrt(2*T)];
end
z = quatmultiply(qacc',qmag');

end

function J=CovEstim(qa,qm,a,l,m)
if iscolumn(l)
    l=l';
end
k=sqrt(1+a(3)); y=l(1:2)*l(1:2)';
b1=sqrt(y+l(1)*sqrt(y)); b2=sqrt(y-l(1)*sqrt(y));
if a(3)>=0
   dqdf1=[qm(1) 0 0 -qm(4) qa(1) -qa(2) qa(3) 0
       0 qm(1) qm(4) 0 qa(2) qa(1) 0 qa(3)
       0 -qm(4) qm(1) 0 qa(3) 0 qa(1) -qa(2)
       qm(4) 0 0 qm(1) 0 -qa(3) qa(2) qa(1)];
   if l(1)>=0
       df1df2=(2*sqrt(2))^-1*...
           [0 0 1/k 0 0 0
           0 -2/k a(2)/k^3 0 0 0
           2/k 0 -a(1)/k^3 0 0 0
           0 0 0 0 0 0
           0 0 0 l(2)^2/(b1*y) prod(l(1:2))/(b1*y) 0
           0 0 0 0 0 0
           0 0 0 0 0 0
           0 0 0 -l(2)*b1/y^(3/2) l(1)*b1/y^(3/2) 0];
   else
       df1df2=(2*sqrt(2))^-1*...
           [0 0 1/k 0 0 0
           0 -2/k a(2)/k^3 0 0 0
           2/k 0 -a(1)/k^3 0 0 0
           0 0 0 0 0 0
           0 0 0 l(2)*b2/y^(3/2) l(1)*b2/y^(3/2) 0
           0 0 0 0 0 0
           0 0 0 0 0 0
           0 0 0 -l(2)^2/(b2*y) prod(l(1:2))/(b2*y) 0];
   end
   df2du=[eye(3) zeros(3)
       m(3)-(2*a(1)*m(1)+a(2)*m(2))/k^2 -a(1)*m(2)/k^2 a(1)*(a(1)*m(1)+a(2)*m(2))/k^4 1-a(1)^2/k^2 -prod(a(1:2))/k^2 a(1)
       -a(2)*m(1)/k^2 m(3)-(a(1)*m(1)+2*a(2)*m(2))/k^2 a(2)*(a(1)*m(1)+a(2)*m(2))/k^4 -prod(a(1:2))/k^2 1-a(2)^2/k^2 a(2)
       -m(1) -m(2) m(3) -a(1) -a(2) a(3)];
    
else
    k=sqrt(1-a(3));
    dqdf1=[qm(1) 0 0 -qm(4) qa(1) -qa(2) 0 -qa(4)
        0 qm(1) qm(4) 0 qa(2) qa(1) -qa(4) 0
        0 -qm(4) qm(1) 0 0 qa(4) qa(1) -qa(2)
        qm(4) 0 0 qm(1) qa(4) 0 qa(2) qa(1)];
    if l(1)>=0
        df1df2=(2*sqrt(2))^-1*...
            [0 -2/k -a(2)/k^3 0 0 0
            0 0 -1/k 0 0 0
            0 0 0 0 0 0
            2/k 0 a(1)/k^3 0 0 0
            0 0 0 l(2)^2/(b1*y) -prod(l(1:2))/(b1*y) 0
            0 0 0 0 0 0
            0 0 0 0 0 0
            0 0 0 -l(2)*b1/y^(3/2) l(1)*b1/y^(3/2) 0];
    else
        df1df2=(2*sqrt(2))^-1*...
            [0 -2/k -a(2)/k^3 0 0 0
            0 0 -1/k 0 0 0
            0 0 0 0 0 0
            2/k 0 a(1)/k^3 0 0 0
            0 0 0 l(2)*b2/y^(3/2) -l(1)*b2/y^(3/2) 0
            0 0 0 0 0 0
            0 0 0 0 0 0
            0 0 0 -l(2)^2/(b2*y) prod(l(1:2))/(b2*y) 0];
    end
    df2du=[eye(3) zeros(3)
        m(3)-(2*a(1)*m(1)-a(2)*m(2))/k^2 a(1)*m(2)/k^2 a(1)*(-a(1)*m(1)+a(2)*m(2))/k^4 1-a(1)^2/k^2 -prod(a(1:2))/k^2 a(1)
        -a(2)*m(1)/k^2 m(3)-(a(1)*m(1)-2*a(2)*m(2))/k^2 a(2)*(-a(1)*m(1)+a(2)*m(2))/k^4 -prod(a(1:2))/k^2 -1+a(2)^2/k^2 a(2)
        m(1) -m(2) -m(3) a(1) -a(2) -a(3)];
end
J=dqdf1*df1df2*df2du;
end