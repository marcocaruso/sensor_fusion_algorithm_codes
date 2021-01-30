function [errAbs1,errAbs2] = EXPLORE_SEL(s1,s2,tauAcc_, tauMag_, zeta, accRating, Qs_, qin1, qin2, indMov, fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(tauAcc_);
M = length(tauMag_);

errAbs1 = zeros(N,M);
errAbs2 = zeros(N,M);

for kk=1:M
    
    tauMag=tauMag_(kk);
    
    for ii=1:N
        
        tauAcc=tauAcc_(ii);
        
        q1 = runOriEstIMU(s1(:,2:4), s1(:,5:7), s1(:,8:10), fs, tauAcc, tauMag, zeta, accRating, qin1);
        q0=quatconj(repmat(q1(1,:),size(q1,1),1));
        sel1=quatmultiply(q0,q1);
        
        q2 = runOriEstIMU(s2(:,2:4), s2(:,5:7), s2(:,8:10), fs, tauAcc, tauMag, zeta, accRating, qin2);
        q0=quatconj(repmat(q2(1,:),size(q2,1),1));
        sel2=quatmultiply(q0,q2);
        
        
        err1 = quatmultiply(quatconj(sel1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(sel2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
        
    end
end
end
