function [R] = quatrotmatr(q)
%this function returns the rotation matrix based on the rotation quaternion
%the scalar part is located at the begin of the quaternion
% $ Version: 1 $
% CODE      by:                 Marco Caruso (Politecnico di Torino,
% Torino, IT) 2017 September 11
% COMMENTS  by:                 Code author             2017 September 11
% OUTPUT    tested by:          Code author             2017 September 11
if size(q,2)~=4
    q=transpose(q);
end

R=zeros(3,3,size(q,1));
for i=1:size(q,1)
    q_rot=q(i,:);
    q0=q_rot(1);
    q1=q_rot(2);
    q2=q_rot(3);
    q3=q_rot(4);
    
    rot_matr=[q0^2+q1^2-q2^2-q3^2, 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);...
        2*(q1*q2+q0*q3), q0^2-q1^2+q2^2-q3^2, 2*(q2*q3-q0*q1);...
        2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), q0^2-q1^2-q2^2+q3^2];
    R(:,:,i)=rot_matr;
end
end

