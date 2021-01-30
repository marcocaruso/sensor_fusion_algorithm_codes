function [errAbs1,errAbs2] = EXPLORE_MAD(s1, s2, beta_, Qs_, qin1, qin2, indMov, fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(beta_);

errAbs1=zeros(N,1);
errAbs2=zeros(N,1);

for ii = 1:length(beta_)
    
    beta=beta_(ii);
    
    q1=zeros(size(s1,1),4);
    q1or = MadgwickAHRS('SamplePeriod', 1/fs);
        
    if beta==0
        q1(1,:)=[1 0 0 0];
        q1or.Quaternion=[1 0 0 0];
    else
        q1(1,:)=qin1;
        q1or.Quaternion=qin1;
    end
    
    q1or.Beta=beta;
    for jj = 2:size(s1,1)
        q1or.Update(s1(jj,5:7), -s1(jj,2:4), s1(jj,8:10));	% gyroscope units must be radians
        q1(jj, :) = q1or.Quaternion;
    end
    
    q0=quatconj(repmat(q1(1,:),size(q1,1),1));
    mad1=quatmultiply(q0,q1);
    
    q2=zeros(size(s1,1),4);
    q2or = MadgwickAHRS('SamplePeriod', 1/fs);
    
    if beta==0
        q2(1,:)=[1 0 0 0];
        q2or.Quaternion=[1 0 0 0];
    else
        q2(1,:)=qin2;
        q2or.Quaternion=qin2;
    end
    
    q2or.Beta=beta;
    for jj = 2:size(s1,1)
        q2or.Update(s2(jj,5:7), -s2(jj,2:4), s2(jj,8:10));	% gyroscope units must be radians
        q2(jj, :) = q2or.Quaternion;
    end
    

    q0=quatconj(repmat(q2(1,:),size(q1,1),1));
    mad2=quatmultiply(q0,q2);
    
    err1 = quatmultiply(quatconj(mad1),Qs_);
    axa1 = quat2axang(err1);
    errAbs1(ii) = rad2deg(rms(axa1(indMov,end)));
    
    err2 = quatmultiply(quatconj(mad2),Qs_);
    axa2 = quat2axang(err2);
    errAbs2(ii) = rad2deg(rms(axa2(indMov,end)));
   
end

end

