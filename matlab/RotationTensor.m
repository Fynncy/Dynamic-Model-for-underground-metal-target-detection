function [Rt] = RotationTensor(theta,phi,psi)
%EULOR ����ŷ���Ǽ���ŷ������
%   ���
%   Rt:ŷ������
%   ���룺
%   theta,������
%   phi,��ת��
%   psi,�����
Rt=zeros(3,3);

% ��̬ת�����󣬲ο�ŷ��������Ǵ�����PPT
% Rt(1,1)=cos(phi)*cos(psi)-sin(phi)*sin(theta)*sin(psi);
% Rt(1,2)=cos(phi)*sin(psi)+sin(phi)*sin(theta)*cos(psi);
% Rt(1,3)=-cos(theta)*sin(phi);
% Rt(2,1)=-cos(theta)*sin(psi);
% Rt(2,2)=cos(theta)*cos(psi);
% Rt(2,3)=sin(theta);
% Rt(3,1)=sin(phi)*cos(psi)+cos(phi)*sin(theta)*sin(psi);
% Rt(3,2)=sin(psi)*sin(phi)-cos(phi)*sin(theta)*cos(psi);
% Rt(3,3)=cos(theta)*cos(phi);

% ��̬ת������2��ŷ�������¸��£�ת��˳��Ϊ�����psi��������theta����ת��phi
% sf=sin(psi);
% cf=cos(psi);
% st=sin(theta);
% ct=cos(theta);
% sg=sin(phi);
% cg=cos(phi);
% Rt(1,1)=cg;
% Rt(1,2)=0;
% Rt(1,3)=-sg;
% Rt(2,1)=st*sg;
% Rt(2,2)=ct;
% Rt(2,3)=st*cg;
% Rt(3,1)=ct*sg;
% Rt(3,2)=-st;
% Rt(3,3)=ct*cg;
% ��̬ת������3���������ҵ�
% th1=theta;
% th2=phi;
% th3=psi;
% R1=[cos(th1),-sin(th1),0;sin(th1),cos(th1),0;0,0,1];
% R2=[1,0,0;0,cos(th2),-sin(th2);0,sin(th2),cos(th2)];
% R3=[cos(th3),-sin(th3),0;sin(th3),cos(th3),0;0,0,1];
% Rt=R1*R2*R3;

% ��̬ת������4����[1] BELL T, COLLINS L. Handheld UXO Sensor Improvements to Facilitate UXO / Clutter Discrimination Volume 1[J]. Work, 2007,  1(November).��¼
R1=[cos(psi),sin(psi),0;-sin(psi),cos(psi),0;0,0,1];   %lw sz 3-28 29 30 31  %lw wz 3-14 15 16
R2=[cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)];
R3=[1,0,0;0,cos(phi),sin(phi);0,-sin(phi),cos(phi)];
Rt=R3*R2*R1;

end

