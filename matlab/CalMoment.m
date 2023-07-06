function [Moment] = CalMoment(MagPolar,Hxyz,theta,phi,psi)
%�������������弫�����Լ�ƫת��
%   ���: 
%       Moment,3x1��������ΪMx,My,Mz
%   ����: 
%       MagPolar������ż�����,betax,betay,betaz
%       Hxyz,3x1��������ΪHx,Hy,Hz
%       theta��������
%       phi����ת��
%       psi������ǣ���ֱ���õ���Գ�������û�еģ�

m_Polar=diag(MagPolar);
Hxyz=Hxyz(:);
Rt=RotationTensor(theta,phi,psi);

Moment=Rt*m_Polar*Rt'*Hxyz;

end

