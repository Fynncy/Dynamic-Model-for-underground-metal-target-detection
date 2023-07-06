function [B] = HAlpha_M(detector,Pk,Alpha,theta,phi)
%������Ȧ����������ģ�Ͳ���������γ� ���������۲ο���֪���������-2.3�ڣ�
%   ���: 
%       B�����γ�
%   ����: 
%       detector������������
%       pk,������λ��
%       Alpha��1*9����������Ϊx,y,z��M11,M22,M33,M12,M13,M23 

  v_r=Alpha(1:3);
  v_M=Alpha(4:9);
  v_r=v_r(:);
  v_M=v_M(:);
  Pk=Pk(:);
  [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3),theta,phi);      %�����ʽ����Gk
  %Bp=FirstField(detector.I,detector.R,detector.F,v_r-Pk);   % ��������Bp %lw
  %ע�͵��˱�����У�
  Bp=FirstField_wz_matrix(detector.I,detector.R,(v_r-Pk)',theta,phi);
  Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...   % ����Wk
      0 Bp(2) 0 Bp(1) 0 Bp(3);...
      0 0 Bp(3) 0 Bp(1) Bp(2)];
   B=Gk*Wk*v_M;

end

