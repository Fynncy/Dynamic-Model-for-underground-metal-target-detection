function [Beta,Angles] = GetBetaAndAngle(v_M)
%�����������ݣ�����̽��Ŀ���ʼ״̬
%   ���: 
%       v_inital����ʼ״̬��[x,y,z,M11,M22,M33,M12,M13,M23] �����������-�����˲����ٵĵ��̽���㷨�о�
%   ����: 
%       v_M��M11,M22,M33,M12,M13,M23
%       
%% ����v_M������ֵ������Ἣ�����Լ���̬��
M=[v_M(1) v_M(4) v_M(5);...
   v_M(4) v_M(2) v_M(6);...
   v_M(5) v_M(6) v_M(3)];
[V,D]=eig(M);                          %lw �����A��ȫ������ֵ�����ɶԽ���D������A��������������V��������
Beta=sort([ D(3,3) D(2,2) D(1,1)]);
Angles(1)=asin(-real(V(1,1)));			%v_M�п����и���������������ֻ��ȡʵ��
Angles(2)=asin(real(V(2,1))/cos(Angles(1)));
Angles(3)=asin(real(V(1,2))/cos(Angles(2)));
Angles=Angles*180/pi;
Angles(Angles<0)=-1*Angles(Angles<0);		%��Ϊ�Ƕ�ֻ����0

end

