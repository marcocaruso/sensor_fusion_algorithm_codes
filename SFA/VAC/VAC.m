function [qvac] = VAC(XX, t, gain_acc_static, gain_mag, th1a, th2a, biasalpha)

% "Keeping a Good Attitude: A Quaternion-Based Orientation Filter for IMUs and MARGs" (Roberto G. Valenti, Ivan Dryanovski and Jizhong Xiao)

% --------------- INPUT ---------------
% XX              = (:,2:4) Accelerometer (m/s^2), (:,5:7) Gyroscope (rad/s), (:,8:10) Magnetometer (a.u.) normalized units
% t               = Nx1 time vector (s)
% gain_acc_static = 1x1 (a.u.) weight to assign to the accelerometer measurements
% thAcc           = 1x1 (a.u.) first threshold value for accelerometer signal
% thMag           = 1x1 (a.u.) second threshold value for accelerometer signal
% biasalpha       = 1x1 (a.u.) weight to remove the gyroscope bias 

% --------------- OUTPUT ---------------
% qvac            = Nx4 [qw qx qy qz]

% Implementation by Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 25/05/2020

dt = mean(diff(t));

g = 9.81;

%Accelerometer
ax = -XX(:,2); ay = -XX(:,3); az = -XX(:,4);

%Gyroscope
wx = XX(:,5); wy = XX(:,6); wz = XX(:,7);

%Magnetometer
hx = XX(:,8); hy = XX(:,9); hz = XX(:,10);

qvac = zeros(length(t),4);

bx = 0; by = 0; bz = 0;

% Thresholds
% gain_acc_static = 0.01;
% gain_mag = 0.01;
% th1a = 0.1; %accelerometer: threshold for adaptive gain (lower limit)
% th2a = 0.99; %accelerometer: threshold for adaptive gain (higher limit)
% biasalpha = 0.01;

kwthr = 0.2;
kaccthr = 0.1;
kdeltawthr = 0.01;
thr_int = 0.9;


for i = 1:length(t)
    
    if i==1
        %first step: initialization
        %ASSUMPTION: THE SYSTEM IS IN A STEADY CONFIGURATION, the gyro
        %output is ignored and the attitude is only determined by signal
        %from accelerometer and magnetometer
        
        n  = norm([ax(i),ay(i),az(i)]);
        axn = ax(i)/n;
        ayn = ay(i)/n;
        azn = az(i)/n;
        
        if azn>=0
            qacc = [sqrt((azn+1)/2) -ayn/sqrt(2*(azn+1)) axn/sqrt(2*(azn+1)) 0]';
        else
            qacc = [-ayn/sqrt(2*(1-azn)) sqrt((1-azn)/2) 0 axn/sqrt(2*(1-azn))]';
        end
        
        hxn= hx(i);
        hyn= hy(i);
        hzn= hz(i);
        
        l = quatrotmatr(qacc)'*[hxn;hyn;hzn];
        
        T = l(1)^2+l(2)^2;
        
        qmag = [sqrt(T+l(1)*sqrt(T))/sqrt(2*T);...
            0;0;...
            l(2)/(sqrt(2)*sqrt(T+l(1)*sqrt(T)))];
        
        q_in = quatmultiply(qacc',qmag');
        
        qvac(1,:) = q_in;
        
    else
        
        %PREDICTION STEP
        %the orientation is determined from gyro output and. Here the
        %quaternion that describe the rate of change of orientation is
        %calculated
        
        steady = true;
        
        if abs(sqrt(ax(i)^2+ay(i)^2+az(i)^2)-g) > kaccthr
            steady=false;
        end
        
        if abs(wx(i)-wx(i-1))> kdeltawthr || abs(wy(i)-wy(i-1))> kdeltawthr || abs(wz(i)-wz(i-1))> kdeltawthr
            steady=false;
        end
        
        if abs(wx(i)-bx)> kwthr || abs(wy(i)-by)> kwthr || abs(wz(i)-bz) > kwthr
            steady=false;
        end
        
        if(steady)
            bx = bx+biasalpha*(wx(i)-bx);
            by = by+biasalpha*(wy(i)-by);
            bz = bz+biasalpha*(wz(i)-bz);
        end
        
        
        wx_un = wx(i)-bx;
        wy_un = wy(i)-by;
        wz_un = wz(i)-bz;
        
        qw = qvac(i-1,:) - 0.5*dt*quatmultiply([0,wx_un,wy_un,wz_un], qvac(i-1,:));
        
        
        qw = quatnormalize(qw);
        
        %CORRECTION STEP
        %1) Accelerometer-Based Correction
        %- determine the predicted gravity and calculate the delta
        %quaternion
        
        n = norm([ax(i),ay(i),az(i)]);
        axn = ax(i)/n;
        ayn = ay(i)/n;
        azn = az(i)/n;
        
        gp = quatrotmatr(qw)'*[axn; ayn; azn];
        
        dqacc = [sqrt((gp(3)+1)/2); -gp(2)/sqrt(2*(gp(3)+1));...
            gp(1)/sqrt(2*(gp(3)+1)); 0];
        
        %Adaptive Gain (accelerometer)
        
        err_s = abs(n-g)/g;
        
        if err_s < th1a
            
            factor = 1;
            
        elseif err_s > th1a && err_s < th2a
            
            factor = 1-((err_s-th1a)/(th2a-th1a));
            
        else
            
            factor = 0;
            
        end
        
        gain_acc = gain_acc_static*factor;
        
        
        if dqacc(1) < thr_int %if it is true compute the Spherical Interpolation
            ang = acos(dqacc(1));
            dqacc = sin((1-gain_acc)*ang)/sin(ang) * [1; 0; 0; 0] + ...
                sin(gain_acc*ang)/sin(ang)*dqacc;
        else
            dqacc = ((1-gain_acc)*[1; 0; 0; 0] + gain_acc*dqacc)./norm(dqacc);
        end
        
        dqacc = quatnormalize(dqacc');
        qa = quatmultiply(qw,dqacc);
        
        
        %2) Magnetometer-based correction
        hxn = hx(i);
        hyn = hy(i);
        hzn = hz(i);
        
        lq = quatrotmatr(qa)'*[hxn;hyn;hzn];
        T = lq(1)^2+lq(2)^2;
        
        dqmag = [sqrt(T+lq(1)*sqrt(T))/sqrt(2*T); 0; 0;...
            lq(2)/sqrt(2*(T+lq(1)*sqrt(T)))];
        
        if dqmag(1)<thr_int %if it is true compute the Spherical Interpolation
            ang = acos(dqmag(1));
            dqmag = sin((1-gain_mag)*ang)/sin(ang) * [1;0;0;0] + ...
                sin(gain_mag*ang)/sin(ang)*dqmag;
        else
            dqmag = ((1-gain_mag)*[1;0;0;0] + gain_mag*dqmag)./norm(dqmag);
        end
        
        qvac(i,:) = quatnormalize(quatmultiply(qa,dqmag'));
        
    end
    
end

end

