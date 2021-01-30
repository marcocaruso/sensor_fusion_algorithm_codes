function [errAbs1,errAbs2] = EXPLORE_LIG(s1, s2, stdGyr_, stdAcc, stdMag, ca, cb_, cm, Qs_, init1, init2, indMov, localMagField, fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(cb_);
M = length(stdGyr_);

errAbs1=zeros(N,M);
errAbs2=zeros(N,M);


Ts = 1/fs;
indQ = 1000;

for kk = 1:M
    
    stdGyro = stdGyr_(kk);
    
    for ii=1:N
        
        cb = cb_(ii);
        
        params_acc = [stdGyro stdAcc [1 1 1].*ca [1 1 1].*cb]';
        params_mag = [stdGyro stdMag [1 1 1].*(cm) [1 1 1].*(cm)]';
        
        
        x_init_acc = [init1.acc_in'; zeros(3, 1)]; % initial conditions from the first nsamp samples
        params_a   = [Ts; init1.g; params_acc; x_init_acc];
        
        x_init_mag = [init1.mag_in'; zeros(3, 1)]; % initial conditions from the first nsamp samples
        params_h   = [Ts; init1.h'; params_mag; x_init_mag];
        clear lkf
        y_acc = lkf_loop(params_a, s1(:,2:4)', s1(:,5:7)');
        clear lkf
        y_mag = lkf_loop(params_h, s1(:,8:10)'*localMagField/init1.normh, s1(:,5:7)');
        
        
        q1 = quat_lkf(y_acc, y_mag, init1.gn, init1.hn');
        q1 = correctQuat(q1);
        
        q0 = quatconj(repmat(q1(indQ,:),size(q1,1),1));
        lig1 = quatmultiply(q0,q1);
        
        
        x_init_acc = [init2.acc_in'; zeros(3, 1)]; % initial conditions from the first nsamp samples
        params_a   = [Ts; init2.g; params_acc; x_init_acc];
        
        x_init_mag = [init2.mag_in'; zeros(3, 1)]; % initial conditions from the first nsamp samples
        params_h   = [Ts; init2.h'; params_mag; x_init_mag];
        
        clear lkf
        y_acc = lkf_loop(params_a, s2(:,2:4)', s2(:,5:7)');
        clear lkf
        y_mag = lkf_loop(params_h, s2(:,8:10)'*localMagField/init2.normh, s2(:,5:7)');
        
        
        q2 = quat_lkf(y_acc, y_mag, init2.gn, init2.hn');
        q2 = correctQuat(q2);
        
        q0 = quatconj(repmat(q2(indQ,:),size(q2,1),1));
        lig2 = quatmultiply(q0,q2);
        
        
        err1 = quatmultiply(quatconj(lig1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(lig2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
        
    end
end
end