function [Hx,Hy,Hz] = SecondField(Coil,Target,RCoil)
%������Ȧ������Ŀ��������۲�ƽ�������������γ�
%   ���: 
%       Hx,Hy,Hz,����Ϊ���γ���x,y,z���ϵķ���
%   ����: 
%       Coil����Ȧ����
%       Target��Ŀ���������
%       Plane���۲��


Mxyz=CalMoment(Target.MagPolar,FirstField(Coil.I,Coil.R,Coil.f,Target.Postion-Coil.Postion),Target.Theta,Target.Phi,Target.Psi);

MM=[Target.Postion(:) Mxyz(:)];   
[Hx,Hy,Hz]=HFieldModel(MM,RCoil.Postion(1),RCoil.Postion(2),RCoil.Postion(3),Target.Theta,Target.Phi);
end

