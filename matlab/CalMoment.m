function [Moment] = CalMoment(MagPolar,Hxyz,theta,phi,psi)
%给定主场和物体极化率以及偏转角
%   输出: 
%       Moment,3x1矩阵，依次为Mx,My,Mz
%   输入: 
%       MagPolar，三轴磁极化率,betax,betay,betaz
%       Hxyz,3x1矩阵，依次为Hx,Hy,Hz
%       theta，俯仰角
%       phi，滚转角
%       psi，航向角（垂直放置的轴对称物体是没有的）

m_Polar=diag(MagPolar);
Hxyz=Hxyz(:);
Rt=RotationTensor(theta,phi,psi);

Moment=Rt*m_Polar*Rt'*Hxyz;

end

