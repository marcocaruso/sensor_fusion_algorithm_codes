function [errAbs1,errAbs2] = EXPLORE_VAK(s1,s2,stdGyr_, stdAcc_, stdMag, Qs_,indMov,fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(stdGyr_);
M = length(stdAcc_);

errAbs1=zeros(N,M);
errAbs2=zeros(N,M);


t = 0:1/fs:(size(s1,1)-1)/fs;


for kk=1:M
    
    stdAcc = stdAcc_(kk);
    
    for ii=1:N
        
        stdGyro = stdGyr_(ii);
        
        q1 = VAK(s1, t, stdAcc, stdGyro, stdMag);
        q1 = quatconj(q1);
        q0 = quatconj(repmat(q1(1,:),size(q1,1),1));
        vak1 = quatmultiply(q0,q1);
        
        q2 = VAK(s2, t, stdAcc, stdGyro, stdMag);
        q2 = quatconj(q2);
        q0 = quatconj(repmat(q2(1,:),size(q2,1),1));
        vak2 = quatmultiply(q0,q2);
        
        
        err1 = quatmultiply(quatconj(vak1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(vak2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
        
    end

end
end
