function [Z,A] = jaccsd_M(detector,Pk,Alpha,theta,phi)
%给定线圈参数，根据模型参数计算二次场在当前位置的雅格比行列式（只关于位置）
%   输出: 
%       Z，当前位置函数值
%       A，雅格比行列式
%   输入: 
%       detector，线圈参数
%       Alpha，1*9向量，依次为x,y,z,M11,M22,M33,M12,M13,M23
%       Pk,探测器位置
%lw     Scene, 添加新参数
Z=HAlpha_M(detector,Pk,Alpha,theta,phi); %lw 添加新参数Scene
gx=zeros(1,length(Alpha));
gy=zeros(1,length(Alpha));
gz=zeros(1,length(Alpha));
Factor=eps^(1/3);
for i=1:length(Alpha)
    if(abs(Alpha(i))<eps)
        DifferenceInterval=Factor;
    else
%         DifferenceInterval=abs(Alpha(i))*Factor;
        DifferenceInterval=Factor;
    end
  
    Alpha(i)=Alpha(i)+DifferenceInterval;       % 第i个参数偏移DifferenceInterval
    [B1]=HAlpha_M(detector,Pk,Alpha,theta,phi);
    Alpha(i)=Alpha(i)-2*DifferenceInterval;     % 第i个参数偏移-DifferenceInterval
    [B2]=HAlpha_M(detector,Pk,Alpha,theta,phi);
    Alpha(i)=Alpha(i)+DifferenceInterval;       % 恢复第i个参数的值
    gx(i)=((B1(1)-B2(1)))/(2*DifferenceInterval);
    gy(i)=((B1(2)-B2(2)))/(2*DifferenceInterval);
    gz(i)=((B1(3)-B2(3)))/(2*DifferenceInterval);
end
A=[gx;gy;gz];
end

