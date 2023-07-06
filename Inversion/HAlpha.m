function [Hx,Hy,Hz] = HAlpha(Coil,Alpha)
%给定线圈参数，根据模型参数计算二次场
%   输出: 
%       Hx,Hy,Hz,依次为二次场在x,y,z轴上的分量
%   输入: 
%       Coil，线圈
%       Alpha，1*9向量，依次为x,y,z,θ，φ，ψ，βx,βy,βz
Postion=Alpha(1:3);
Mxyz=CalMoment(Alpha(7:9),FirstField(Coil.I,Coil.R,Coil.f,Postion(:)-Coil.Postion(:)),Alpha(4),Alpha(5),Alpha(6));

MM=[Postion(:) Mxyz(:)];   
[Hx,Hy,Hz]=HFieldModel(MM,Coil.Postion(1),Coil.Postion(2),Coil.Postion(3),Alpha(4),Alpha(5));

end

