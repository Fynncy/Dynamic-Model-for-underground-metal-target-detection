function [Hx,Hy,Hz] = SecondField(Coil,Target,RCoil)
%给定线圈参数，目标参数，观测平面参数，计算二次场
%   输出: 
%       Hx,Hy,Hz,依次为二次场在x,y,z轴上的分量
%   输入: 
%       Coil，线圈参数
%       Target，目标物体参数
%       Plane，观测点


Mxyz=CalMoment(Target.MagPolar,FirstField(Coil.I,Coil.R,Coil.f,Target.Postion-Coil.Postion),Target.Theta,Target.Phi,Target.Psi);

MM=[Target.Postion(:) Mxyz(:)];   
[Hx,Hy,Hz]=HFieldModel(MM,RCoil.Postion(1),RCoil.Postion(2),RCoil.Postion(3),Target.Theta,Target.Phi);
end

