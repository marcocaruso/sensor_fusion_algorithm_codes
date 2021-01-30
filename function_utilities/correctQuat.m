function [Q] = correctQuat(Q)
% Qc=zeros(size(Q));
% Qc(1,:)=Q(1,:);
for i=2:size(Q,1)
    ax1=quat2axang(Q(i-1,:));
    ax2=quat2axang(Q(i,:));
    an=real(acosd(dot(ax1(1:3),ax2(1:3))));
    copia(i)=an;
    if an<150
        Q(i,:)=Q(i,:);
    else
        Q(i,:)=-Q(i,:);
    end
end
end

