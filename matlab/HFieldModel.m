function [Bx,By,Bz,G] = HFieldModel(MM,x,y,z,theta,phi)
%HFIELD �����ż������x,y,z���ĴŸ�Ӧǿ��B����λΪnT
%   �����
%       Bx,By,Bz���ų���x,y,z��������ķ���
%       G�����ֺ���3*3����
%   ���룺
%       MM,�žأ�3*2���󣬵�һ��Ϊλ�ò������ڶ���Ϊx,y,zλ�ã�
%       x,y,z,�۲��λ��
%  
T=zeros(3,3);
m_p=MM(:,1);    % λ��
%m=10^-7*MM(:,2);      % ż����
m=MM(:,2);         % ż����
G=zeros(3);            % ���ֺ���

rt = RotationTensor(theta,phi,0);   %lw ����ȱ�Σ���������
x=x-m_p(1);
y=y-m_p(2);
z=z-m_p(3);

r=sqrt(x.^2+y.^2+z.^2);
r5=1./r.^5;
r3=1./r.^3;

G(1,1)=3*x.*x.*r5-r3;
G(1,2)=3*x.*y.*r5;
G(1,3)=3*x.*z.*r5;
G(2,1)=3*x.*y.*r5;
G(2,2)=3*y.*y.*r5-r3;
G(2,3)=3*y.*z.*r5;
G(3,1)=3*x.*z.*r5;
G(3,2)=3*y.*z.*r5;
G(3,3)=3*z.*z.*r5-r3;
G=100*G;           % Ϊ��֤���յ�λΪnT
Bx=m(1)*G(1,1)+m(2)*G(1,2)+m(3)*G(1,3);  %lw wz 3-5
By=m(1)*G(2,1)+m(2)*G(2,2)+G(2,3)*m(3);
Bz=G(3,1)*m(1)+G(3,2)*m(2)+G(3,3)*m(3);
% Bx=abs(Bx);
% By=abs(By);
% Bz=abs(Bz);

B = rt' * [Bx,By,Bz]';   %lw wz 3-21
% B = abs(B);
Bx = B(1);
By = B(2);
Bz = B(3);

end

