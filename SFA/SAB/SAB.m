function qsab = SAB(acc, gyr, mag, t, stdAcc, stdGyro, stdMag, thAcc, thMag, stdBiasMag_in, qin)
% Sabatini 2011 - Estimating three dimensional orientation of human body parts by inertial-magnetic sensing (Sensors, 2011)

% --------------- INPUT ---------------
% acc            = Nx3 (m/s^2)
% gyr            = Nx3 (rad/s)
% mag            = Nx3 (a.u.) normalized units
% t              = Nx1 (s)
% stdAcc         = 1x1 (m/s^2) inverse weight to assign to the accelerometer measurements
% stdGyro        = 1x1 (rad/s) inverse weight to assign to the gyroscope measurements
% stdMag         = 1x1 (a.u.)  inverse weight to assign to the magnetometer measurements
% thAcc          = 1x1 (m/s^2) threshold value for accelerometer signal
% thMag          = 1x1 (a.u.) threshold value for magnetometer signal
% stdBiasMag_in  = 1x1 (a.u.) initial value for magnetometer bias covariance matrix
% qin            = 1x4  initial quaternion

% --------------- OUTPUT ---------------
% qsab           = Nx4 [qx qy qz qw], the scalar part is at the END of the quaternion

% Implementation by Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 21/11/2019


% Accelerometer data
ax = -acc(:,1); ay = -acc(:,2); az = -acc(:,3);

% Gyroscope data
wx = gyr(:,1); wy = gyr(:,2); wz = gyr(:,3);

% Magnetometer data
hx = mag(:,1); hy = mag(:,2); hz = mag(:,3);


xPriori=zeros(7,length(t));
xPost=zeros(7,length(t));

[~,~,~,L] = initQuaternion(-acc(1,:),mag(1,:));

% if ~exist('qin','var')
%     [qin,~,~,L] = initialEKFquat(-acc(1,:),mag(1,:));
% else
%     [~,~,~,L] = initialEKFquat(-acc(1,:),mag(1,:));
% end

if isrow(qin)
    qin=transpose(qin);
end

qin=[-qin(2:end); qin(1)]; %da loc a glob
xPost(1:4,1)=qin;

%PROCESS NOISE
SIGMA_g=stdGyro^2*eye(3);

%Constants
g=9.81;
h=[sqrt(L(1).^2+L(2).^2);0;L(3)]; %Earth's magnetic field (global Frame)

dt=mean(diff(t)); %only to initialize the SIGMA_m matrix
SIGMA_m=dt*stdBiasMag_in^2*eye(3);

CSI=[[0 -xPost(3,1) xPost(2,1);
    xPost(3,1) 0 -xPost(1,1)
    -xPost(2,1) xPost(1,1) 0]+xPost(4,1)*eye(3);...
    -xPost(1:3,1)'];

Q=[(dt/2)^2*CSI*SIGMA_g*CSI' zeros(4,3)
    zeros(3,4) SIGMA_m];

Ppost=Q; %posterior initial guess covariance error matrix


warning off
for i=1:length(t)-1
    dt = t(i+1) - t(i);
    %PREDICTION STEP
    omega=0.5*[0 wz(i) -wy(i) wx(i)
        -wz(i) 0 wx(i) wy(i)
        wy(i) -wx(i) 0 wz(i)
        -wx(i) -wy(i) -wz(i) 0];%skew symmetric
    
    F=[eye(4)+omega*dt+0.5*(omega*dt)^2 zeros(4,3)
        zeros(3,4) eye(3)]; %linearized to increase computational speed
    
%     F=[expm(omega*dt) zeros(4,3)
%         zeros(3,4) eye(3)];
    
    %Project the state ahead
    xPriori(:,i)=F*xPost(:,i);
    
    CSI=[[0 -xPost(3,i) xPost(2,i);
        xPost(3,i) 0 -xPost(1,i)
        -xPost(2,i) xPost(1,i) 0]+xPost(4,i)*eye(3);...
        -xPost(1:3,i)']; 
        
    Q=[(dt/2)^2*CSI*SIGMA_g*CSI' zeros(4,3)
        zeros(3,4) SIGMA_m];
    
    %Compute the a priori covariance matrix
    Ppriori= F*Ppost*F'+Q;
    
    %UPDATE STEP
    %Linearize the measurement equation: Jacobian
    q1=xPriori(1,i);
    q2=xPriori(2,i);
    q3=xPriori(3,i);
    q4=xPriori(4,i);
    
    H=[[2*g*[q3 -q4 q1 -q2
        q4 q3 q2 q1
        -q1 -q2 q3 q4];
        2*[q1*h(1)+q3*h(3) -q2*h(1)-q4*h(3) -q3*h(1)+q1*h(3) q4*h(1)-q2*h(3)
        q2*h(1)+q4*h(3) q1*h(1)+q3*h(3) -q4*h(1)+q2*h(3) -q3*h(1)+q1*h(3)
        q3*h(1)-q1*h(3) q4*h(1)-q2*h(3) q1*h(1)+q3*h(3) q2*h(1)+q4*h(3)]
        ] [zeros(3,3);eye(3,3)]];
    
    C=[q1^2-q2^2-q3^2+q4^2 2*(q1*q2+q3*q4) 2*(q1*q3-q2*q4)
        2*(q1*q2-q3*q4) -q1^2+q2^2-q3^2+q4^2 2*(q2*q3+q1*q4)
        2*(q1*q3+q2*q4) 2*(q2*q3-q1*q4) -q1^2-q2^2+q3^2+q4^2];
    
    %Adapt the measurement covariance matrix
    %Accelerometer
    if norm([ax(i);ay(i);az(i)]-C*[0;0;g])<thAcc
        std_acc=stdAcc;
    else
        std_acc=1e100;
    end
    %Magnetometer
    if norm([hx(i);hy(i);hz(i)]-xPriori(5:end,i)-C*h)<thMag
        std_mag=stdMag;
    else
        std_mag=1e100;
    end
    
    
    % Measurement covariance
    R=[std_acc^2*eye(3) zeros(3)
        zeros(3) std_mag^2*eye(3)];
    
    %Compute the Kalman Gain
    K=Ppriori*H'*(H*Ppriori*H'+R)^-1;
    
    %Compute the a posteriori state estimate
    z=[ax(i);ay(i);az(i);hx(i);hy(i);hz(i)]; %measurement vector
    z_predict=[C zeros(3); zeros(3) C]*[0;0;g;h]+[0;0;0;xPriori(5:end,i)];
    xPost(:,i+1)=xPriori(:,i)+K*(z-z_predict);
    
    %Normalize the quaternion
    xPost(1:4,i+1)=xPost(1:4,i+1)/norm(xPost(1:4,i+1));
    
    %Compute the a posteriori covariance matrix
    Ppost=Ppriori-K*H*Ppriori;
end

qsab=xPost(1:4,:)';

end


function [z,qacc,qmag,l] = initialEKFquat(acc,mag)
% Valenti 2015 - Keeping a Good Attitude: A Quaternion-Based Orientation
% Filter for IMUs and MARGs (Sensors, 2015)

% Implementation: Marco Caruso (Politecnico di Torino) 
% Date: 21/11/2019

% --------------- INPUT ---------------
% acc            = 1x3 (m/s^2)
% mag            = 1x3 (a.u.) normalized units

% --------------- OUTPUT ---------------
% z              = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% qacc           = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% qmag           = 1x4 [qw qx qy qz], the scalar part is at the BEGINNING of the quaternion
% l              = 1x3 measured field rotated in the global horizontal plane

n=norm(acc);
ax=acc(1)/n; ay=acc(2)/n; az=acc(3)/n;
if az>=0
    qacc=[sqrt((az+1)/2) -ay/sqrt(2*(az+1)) ax/sqrt(2*(az+1)) 0]';
else
    qacc=[-ay/sqrt(2*(1-az)) sqrt((1-az)/2) 0 ax/sqrt(2*(1-az))]';
end

hx=mag(1); hy=mag(2); hz=mag(3);
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