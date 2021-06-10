function [quaternion] = GUO(XX, stdGyro, stdAcc, stdMag, qin)
%Fast Kalman Filter for Attitude Estimation (S. Guo, J. Wu, Z. Wang, and J.Qian, “Novel MARG-Sensor Orientation Estimation Algorithm Using Fast Kalman Filter,” J. Sensors, vol. 2017, 2017, doi: 10.1155/2017/8542153)

% --------------- INPUT ---------------
% XX             = (:,1) time vector (seconds)
% XX             = (:,2:4) Accelerometer (m/s^2), (:,5:7) Gyroscope (rad/s), (:,8:10) Magnetometer (a.u.) normalized units
% stdAcc         = 1x1 (m/s^2) inverse weight to assign to the accelerometer measurements
% stdGyro        = 1x1 (rad/s) inverse weight to assign to the gyroscope measurements
% stdMag         = 1x1 (a.u.)  inverse weight to assign to the magnetometer measurements
% qin            = 1x4  initial quaternion


% --------------- OUTPUT ---------------
% quaternion     = Nx4 [qw qx qy qz]

%author: Jin Wu
%e-mail: jin_wu_uestc@hotmail.com
% Adapted by Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 30/08/2020



Accelerometer = -XX(:, 2:4);
Gyroscope = XX(:, 5:7);
Magnetometer = XX(:, 8:10);

Pk = 0.001 * eye(4);

if exist('qin','var')
    if ~iscolumn(qin)
        q = transpose(qin);
    else
        q = qin;
    end
else
    q = [1;0;0;0];
end

dt = mean(diff(XX(:,1)));
len=length(Accelerometer(:,1));
quaternion=zeros(len,4);

Sigma_a = stdAcc * eye(3);
Sigma_g = stdGyro * eye(3);
Sigma_m = stdMag * eye(3);

for i=1:len
    q0=q(1);
    q1=q(2);
    q2=q(3);
    q3=q(4);
    
    wx=Gyroscope(i,1);
    wy=Gyroscope(i,2);
    wz=Gyroscope(i,3);
    
    Accelerometer(i,:)=Accelerometer(i,:)./norm(Accelerometer(i,:));
    Magnetometer(i,:)=Magnetometer(i,:)./norm(Magnetometer(i,:));
    
    mD=dot(Accelerometer(i,:),Magnetometer(i,:));
    mN=sqrt(1-mD^2);
    
    omega4=[0,-wx,-wy,-wz;
        wx,0,wz,-wy;
        wy,-wz,0,wx;
        wz,wy,-wx,0];
    
    Phi=eye(4)+dt/2*omega4;
    
    Dk=[q1 q2 q3;
        -q0 -q3 -q2;
        q2 -q0 -q1;
        -q2 q1 -q0];
    Xi=dt*dt/4*Dk*Sigma_g*Dk';
    
    [qy, Jacob]=measurement_quaternion_acc_mag(Accelerometer(i,:),Magnetometer(i,:),[mN,0,mD], q);
    qy=qy./norm(qy);
    
    Eps=Jacob*[Sigma_a,zeros(3,3);zeros(3,3),Sigma_m]*Jacob';
    
    q_=q;
    Pk_ = Pk;
    [q , Pk] = kalman_update(q_, qy, Pk_, Phi, Xi, Eps);
    
    q=q./norm(q);
    quaternion(i,:)=q';
end

end

function [xk , Pk ] = kalman_update(xk_1,yk,Pk_1,Phi_k,Xi_k,Eps_k )
x_=Phi_k*xk_1;
Pk_=Phi_k*Pk_1*Phi_k'+Xi_k;
Gk=Pk_*(inv(Pk_+Eps_k));
Pk=(eye(4)-Gk)*Pk_;
xk=x_+Gk*(yk-x_);
end

function [q,Jacob]=measurement_quaternion_acc_mag(acc,mag,mag_r,q_)

ax=acc(1);  ay=acc(2);  az=acc(3);
mx=mag(1);  my=mag(2);  mz=mag(3);
mN=mag_r(1);            mD=mag_r(3);

q0=q_(1);   q1=q_(2);   q2=q_(3);   q3=q_(4);

q=zeros(4,1);
Jacob=zeros(4,6);

q(1)= (ay*mD*my + (1 + az)*(1 + mN*mx + mD*mz) + ax*(mD*mx - mN*mz))*q0 + ((mD + az*mD - ax*mN)*my + ay*(1 + mN*mx - mD*mz))*q1 + ...
    (ay*mN*my + ax*(-1 + mN*mx + mD*mz) + (1 + az)*(-(mD*mx) + mN*mz))*q2 + (-((ax*mD + mN + az*mN)*my) + ay*(mD*mx + mN*mz))*q3;

q(2)= ((mD - az*mD - ax*mN)*my + ay*(1 + mN*mx + mD*mz))*q0 + (ay*mD*my - (-1 + az)*(1 + mN*mx - mD*mz) + ax*(mD*mx + mN*mz))*q1 + ...
    ((ax*mD + mN - az*mN)*my + ay*(-(mD*mx) + mN*mz))*q2 + (-(ay*mN*my) + ax*(1 - mN*mx + mD*mz) - (-1 + az)*(mD*mx + mN*mz))*q3;

q(3)= (-(ay*mN*my) - ax*(1 + mN*mx + mD*mz) + (-1 + az)*(mD*mx - mN*mz))*q0 + ((-(ax*mD) + mN - az*mN)*my + ay*(mD*mx + mN*mz))*q1 + ...
    (ay*mD*my + (-1 + az)*(-1 + mN*mx + mD*mz) + ax*(mD*mx - mN*mz))*q2 + ((mD - az*mD + ax*mN)*my + ay*(1 - mN*mx + mD*mz))*q3;

q(4)= ax*(q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3)) + (1 + az)*(mD*mx*q1 + mD*my*q2 + q3 + mD*mz*q3 - mN*(my*q0 - mz*q1 + mx*q3)) + ...
    ay*(mN*mz*q0 + mN*my*q1 + q2 - mN*mx*q2 - mD*(mx*q0 + mz*q2 - my*q3));

Jacob(1,1)= -q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
Jacob(1,2)= q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
Jacob(1,3)= q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
Jacob(1,4)= (ax*mD + mN + az*mN)*q0 + ay*mN*q1 + (-((1 + az)*mD) + ax*mN)*q2 + ay*mD*q3;
Jacob(1,5)= ay*mD*q0 + (mD + az*mD - ax*mN)*q1 + ay*mN*q2 - (ax*mD + mN + az*mN)*q3;
Jacob(1,6)= mD*(q0 + az*q0 - ay*q1 + ax*q2) + mN*(-(ax*q0) + q2 + az*q2 + ay*q3);

Jacob(2,1)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(2,2)= q0 + mN*mx*q0 + mD*mz*q0 + mD*my*q1 - mD*mx*q2 + mN*mz*q2 - mN*my*q3;
Jacob(2,3)= -((1 + mN*mx)*q1) - mD*(my*q0 - mz*q1 + mx*q3) - mN*(my*q2 + mz*q3);
Jacob(2,4)= ay*(mN*q0 - mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(mD*q1 - mN*q3);
Jacob(2,5)= mD*(q0 - az*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + (-1 + az)*q2 + ay*q3);
Jacob(2,6)= ay*(mD*q0 + mN*q2) + mD*((-1 + az)*q1 + ax*q3) + mN*(ax*q1 + q3 - az*q3);

Jacob(3,1)= -((1 + mN*mx + mD*mz)*q0) - mD*my*q1 + mD*mx*q2 - mN*mz*q2 + mN*my*q3;
Jacob(3,2)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(3,3)= -q2 - mN*(mz*q0 + my*q1 - mx*q2) + mD*(mx*q0 + mz*q2 - my*q3);
Jacob(3,4)= mD*((-1 + az)*q0 + ay*q1 + ax*q2) - mN*(ax*q0 + q2 - az*q2 + ay*q3);
Jacob(3,5)= ay*(-(mN*q0) + mD*q2) - (-1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);
Jacob(3,6)= mN*(q0 - az*q0 + ay*q1) - ax*(mD*q0 + mN*q2) + mD*((-1 + az)*q2 + ay*q3);

Jacob(4,1)= q1 + mN*mx*q1 + mN*my*q2 + mN*mz*q3 + mD*(my*q0 - mz*q1 + mx*q3);
Jacob(4,2)= q2 + mN*(mz*q0 + my*q1 - mx*q2) - mD*(mx*q0 + mz*q2 - my*q3);
Jacob(4,3)= q3 - mN*(my*q0 - mz*q1 + mx*q3) + mD*(mx*q1 + my*q2 + mz*q3);
Jacob(4,4)= -(ay*(mD*q0 + mN*q2)) + ax*(mN*q1 + mD*q3) + (1 + az)*(mD*q1 - mN*q3);
Jacob(4,5)= (1 + az)*(-(mN*q0) + mD*q2) + ax*(mD*q0 + mN*q2) + ay*(mN*q1 + mD*q3);
Jacob(4,6)= ay*(mN*q0 - mD*q2) + (1 + az)*(mN*q1 + mD*q3) + ax*(-(mD*q1) + mN*q3);

Jacob=Jacob*0.25;

end

