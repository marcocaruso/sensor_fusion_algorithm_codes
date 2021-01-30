function [errAbs1,errAbs2] = EXPLORE_SAB(s1,s2,stdGyr_, stdAcc, stdMag, thAcc_, thMag, stdBiasMag,Qs_,qin1,qin2,indMov,fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(stdGyr_);
M = length(thAcc_);

errAbs1 = zeros(N,M);
errAbs2 = zeros(N,M);


t = 0:1/fs:(size(s1,1)-1)/fs;

for kk = 1:M
    
    thAcc = thAcc_(kk);
 
    for ii=1:N
        
        stdGyro = stdGyr_(ii);
        
        q1 = SAB(s1(:,2:4), s1(:,5:7), s1(:,8:10), t, stdAcc, stdGyro, stdMag, thAcc, thMag, stdBiasMag, qin1);
        q1 = [q1(:,4) q1(:,1:3)];
        q0 = quatconj(repmat(q1(1,:),size(q1,1),1));
        sab1 = quatmultiply(q0,q1);
        
        q2 = SAB(s2(:,2:4), s2(:,5:7), s2(:,8:10), t, stdAcc, stdGyro, stdMag, thAcc, thMag, stdBiasMag, qin2);
        q2 = [q2(:,4) q2(:,1:3)];
        q0 = quatconj(repmat(q2(1,:),size(q2,1),1));
        sab2 = quatmultiply(q0,q2);
        
        
        
        err1 = quatmultiply(quatconj(sab1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(sab2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
     
    end
    
end
end