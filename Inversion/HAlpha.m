function [Hx,Hy,Hz] = HAlpha(Coil,Alpha)
%������Ȧ����������ģ�Ͳ���������γ�
%   ���: 
%       Hx,Hy,Hz,����Ϊ���γ���x,y,z���ϵķ���
%   ����: 
%       Coil����Ȧ
%       Alpha��1*9����������Ϊx,y,z,�ȣ��գ��ף���x,��y,��z
Postion=Alpha(1:3);
Mxyz=CalMoment(Alpha(7:9),FirstField(Coil.I,Coil.R,Coil.f,Postion(:)-Coil.Postion(:)),Alpha(4),Alpha(5),Alpha(6));

MM=[Postion(:) Mxyz(:)];   
[Hx,Hy,Hz]=HFieldModel(MM,Coil.Postion(1),Coil.Postion(2),Coil.Postion(3),Alpha(4),Alpha(5));

end

