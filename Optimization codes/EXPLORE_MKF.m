function [errAbs1,errAbs2] = EXPLORE_MKF(s1,s2,varGyr_, Qs_,qin1,qin2,indMov,fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(varGyr_);

errAbs1=zeros(N,1);
errAbs2=zeros(N,1);

for ii = 1:N
    
    varGyr = varGyr_(ii);
    
    ahrs1 = ahrsfilter;
    ahrs1.SampleRate = fs;
    ahrs1.GyroscopeNoise = varGyr;
    ahrs1.ExpectedMagneticFieldStrength = 45.76;
    
    q1=correctQuat(compact(ahrs1(s1(:,2:4),s1(:,5:7),s1(:,8:10))));
    q0 = quatconj(repmat(q1(1,:),size(q1,1),1));
    mkf1 = quatmultiply(q0,q1);
    
    ahrs2 = ahrsfilter;
    ahrs2.SampleRate = fs;
    ahrs2.GyroscopeNoise = varGyr;
    ahrs2.ExpectedMagneticFieldStrength = 45.76;
    
    q2=correctQuat(compact(ahrs2(s2(:,2:4),s2(:,5:7),s2(:,8:10))));
    q0 = quatconj(repmat(q2(1,:),size(q2,1),1));
    mkf2 = quatmultiply(q0,q2);
    
    
    err1 = quatmultiply(quatconj(mkf1),Qs_);
    axa1 = quat2axang(err1);
    errAbs1(ii) = rad2deg(rms(axa1(indMov,end)));
    
    err2 = quatmultiply(quatconj(mkf2),Qs_);
    axa2 = quat2axang(err2);
    errAbs2(ii) = rad2deg(rms(axa2(indMov,end)));

end
end