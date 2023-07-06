function [Z,A] = jaccsd(Coil,Alpha,MagPolar)
%������Ȧ����������ģ�Ͳ���������γ��ڵ�ǰλ�õ��Ÿ������ʽ
%   ���: 
%       Z����ǰλ�ú���ֵ
%       A���Ÿ������ʽ
%   ����: 
%       Coil����Ȧ
%       Alpha��1*6����������Ϊx,y,z,�ȣ��գ��ף���x,��y,��z
%       MagPolar,1*3����������Ϊ��x,��y,��z
Z=zeros(3,1);
[Z(1),Z(2),Z(3)]=HAlpha(Coil,[Alpha MagPolar]);
gx=zeros(size(Alpha));
gy=zeros(size(Alpha));
gz=zeros(size(Alpha));
Factor=eps^(1/3);
for i=1:length(Alpha)
    if(abs(Alpha(i))<eps)
        DifferenceInterval=Factor;
    else
%         DifferenceInterval=abs(Alpha(i))*Factor;
        DifferenceInterval=Factor;
    end
  
    Alpha(i)=Alpha(i)+DifferenceInterval;       % ��i������ƫ��DifferenceInterval
    [hx1,hy1,hz1]=HAlpha(Coil,[Alpha MagPolar]);
    Alpha(i)=Alpha(i)-2*DifferenceInterval;     % ��i������ƫ��-DifferenceInterval
    [hx2,hy2,hz2]=HAlpha(Coil,[Alpha MagPolar]);
    Alpha(i)=Alpha(i)+DifferenceInterval;       % �ָ���i��������ֵ
    gx(i)=((hx1-hx2))/(2*DifferenceInterval);
    gy(i)=((hy1-hy2))/(2*DifferenceInterval);
    gz(i)=((hz1-hz2))/(2*DifferenceInterval);
end
A=[gx;gy;gz];
end

