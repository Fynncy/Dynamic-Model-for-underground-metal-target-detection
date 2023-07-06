function [Z,A] = jaccsd(Coil,Alpha,MagPolar)
%给定线圈参数，根据模型参数计算二次场在当前位置的雅格比行列式
%   输出: 
%       Z，当前位置函数值
%       A，雅格比行列式
%   输入: 
%       Coil，线圈
%       Alpha，1*6向量，依次为x,y,z,θ，φ，ψ，βx,βy,βz
%       MagPolar,1*3向量，依次为βx,βy,βz
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
  
    Alpha(i)=Alpha(i)+DifferenceInterval;       % 第i个参数偏移DifferenceInterval
    [hx1,hy1,hz1]=HAlpha(Coil,[Alpha MagPolar]);
    Alpha(i)=Alpha(i)-2*DifferenceInterval;     % 第i个参数偏移-DifferenceInterval
    [hx2,hy2,hz2]=HAlpha(Coil,[Alpha MagPolar]);
    Alpha(i)=Alpha(i)+DifferenceInterval;       % 恢复第i个参数的值
    gx(i)=((hx1-hx2))/(2*DifferenceInterval);
    gy(i)=((hy1-hy2))/(2*DifferenceInterval);
    gz(i)=((hz1-hz2))/(2*DifferenceInterval);
end
A=[gx;gy;gz];
end

