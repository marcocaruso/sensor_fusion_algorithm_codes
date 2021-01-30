function [errAbs1,errAbs2] = EXPLORE_MAH(s1,s2,kp_,ki_,Qs_,qin1,qin2,indMov,fs)

% Loop to compute absolute orientation erro for each unit
% Please refer to the specific SFA for a detailed explanation of the
% variables
% Authors: Marco Caruso (marco.caruso@polito.it)
% PolitoBIOMed Lab – Biomedical Engineering Lab and Department of Electronics and Telecommunications, Politecnico di Torino, Torino, Italy; 
% Last modified: 31/01/2020

N = length(kp_);
M = length(ki_);

errAbs1=zeros(N,M);
errAbs2=zeros(N,M);


for kk = 1:M
    ki=ki_(kk);
    for ii=1:N
        
        kp=kp_(ii);
        
        q1=zeros(size(s1,1),4);
        q1or = MahonyAHRS('SamplePeriod', 1/fs);
        
        if kp==0
            q1(1,:)=[1 0 0 0];
            q1or.Quaternion=[1 0 0 0];
        else
            q1(1,:)=qin1;
            q1or.Quaternion=qin1;
        end
        
        q1or.Kp=kp;
        q1or.Ki=ki;
        
        for jj = 2:size(s1,1)
            q1or.Update(s1(jj,5:7), -s1(jj,2:4), s1(jj,8:10));	% gyroscope units must be radians
            q1(jj, :) = q1or.Quaternion;
        end
        
        q0=quatconj(repmat(q1(1,:),size(q1,1),1));
        mah1=quatmultiply(q0,q1);
        
        q2=zeros(size(s1,1),4);
        q2or = MahonyAHRS('SamplePeriod', 1/fs);
        
        if kp==0
            q2(1,:)=[1 0 0 0];
            q2or.Quaternion=[1 0 0 0];
        else
            q2(1,:)=qin2;
            q2or.Quaternion=qin2;
        end
        
        q2or.Kp=kp;
        q2or.Ki=ki;
        for jj = 2:size(s1,1)
            q2or.Update(s2(jj,5:7), -s2(jj,2:4), s2(jj,8:10));	% gyroscope units must be radians
            q2(jj, :) = q2or.Quaternion;
        end

        q0=quatconj(repmat(q2(1,:),size(q1,1),1));
        mah2=quatmultiply(q0,q2);
        
        
        err1 = quatmultiply(quatconj(mah1),Qs_);
        axa1 = quat2axang(err1);
        errAbs1(ii,kk) = rad2deg(rms(axa1(indMov,end)));
        
        err2 = quatmultiply(quatconj(mah2),Qs_);
        axa2 = quat2axang(err2);
        errAbs2(ii,kk) = rad2deg(rms(axa2(indMov,end)));
        
        
    end
end

end

