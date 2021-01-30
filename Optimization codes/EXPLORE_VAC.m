function [errAbs1,errAbs2] = EXPLORE_VAC(s1,s2, gainAcc_, gainMag_, th1a_, th2a_, biasalpha, Qs_, indMov, fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(gainMag_);
M = length(th2a_);

errAbs1=zeros(N,M);
errAbs2=zeros(N,M);


t = 0:1/fs:(size(s1,1)-1)/fs;

gain_acc_static = gainAcc_;

for kk = 1:M
    
    th2a = th2a_(kk);
    th1a = th1a_;
    
    for ii=1:N
        
        gain_mag = gainMag_(ii);
        
        
        q1 = VAC(s1, t, gain_acc_static, gain_mag, th1a, th2a, biasalpha);
        q1 = quatconj(q1);
        q0 = quatconj(repmat(q1(1,:),size(q1,1),1));
        vac1 = quatmultiply(q0,q1);
        
        q2 = VAC(s2, t, gain_acc_static, gain_mag, th1a, th2a, biasalpha);
        q2 = quatconj(q2);
        q0 = quatconj(repmat(q2(1,:),size(q2,1),1));
        vac2 = quatmultiply(q0,q2);
        
        
        err1 = quatmultiply(quatconj(vac1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(vac2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
        
    end
end
end
