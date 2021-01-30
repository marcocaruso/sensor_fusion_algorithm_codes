function [errAbs1,errAbs2] = EXPLORE_GUO(s1, s2, stdGyr_, stdAcc, stdMag, Qs_, indMov, qin1, qin2, fs)
% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(stdGyr_);
errAbs1 = zeros(N,1);
errAbs2 = zeros(N,1);


for ii=1:N
    
    stdGyro = stdGyr_(ii);
    
    q1 = GUO(s1, stdGyro, stdAcc, stdMag, qin1);
    q0 = quatconj(repmat(q1(5000,:),size(q1,1),1));
    guo1 = quatmultiply(q0,q1);
    
    q2 = GUO(s2, stdGyro, stdAcc, stdMag, qin2);
    q0 = quatconj(repmat(q2(5000,:),size(q2,1),1));
    guo2 = quatmultiply(q0,q2);
    
    
    err1 = quatmultiply(quatconj(guo1),Qs_);
    axa1 = quat2axang(err1);
    errAbs1(ii) = rad2deg(rms(axa1(indMov,end)));
    
    err2 = quatmultiply(quatconj(guo2),Qs_);
    axa2 = quat2axang(err2);
    errAbs2(ii) = rad2deg(rms(axa2(indMov,end)));
    
end
warning on
end
