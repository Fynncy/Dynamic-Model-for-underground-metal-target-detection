function [Z,A] = jaccsd_M(detector,Pk,Alpha,theta,phi)
%������Ȧ����������ģ�Ͳ���������γ��ڵ�ǰλ�õ��Ÿ������ʽ��ֻ����λ�ã�
%   ���: 
%       Z����ǰλ�ú���ֵ
%       A���Ÿ������ʽ
%   ����: 
%       detector����Ȧ����
%       Alpha��1*9����������Ϊx,y,z,M11,M22,M33,M12,M13,M23
%       Pk,̽����λ��
%lw     Scene, ����²���
Z=HAlpha_M(detector,Pk,Alpha,theta,phi); %lw ����²���Scene
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
  
    Alpha(i)=Alpha(i)+DifferenceInterval;       % ��i������ƫ��DifferenceInterval
    [B1]=HAlpha_M(detector,Pk,Alpha,theta,phi);
    Alpha(i)=Alpha(i)-2*DifferenceInterval;     % ��i������ƫ��-DifferenceInterval
    [B2]=HAlpha_M(detector,Pk,Alpha,theta,phi);
    Alpha(i)=Alpha(i)+DifferenceInterval;       % �ָ���i��������ֵ
    gx(i)=((B1(1)-B2(1)))/(2*DifferenceInterval);
    gy(i)=((B1(2)-B2(2)))/(2*DifferenceInterval);
    gz(i)=((B1(3)-B2(3)))/(2*DifferenceInterval);
end
A=[gx;gy;gz];
end

