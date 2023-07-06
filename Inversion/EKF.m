function [m_Alpha] = EKF(Scene,Alpha0)
%�����������ݺ�̽��Ŀ���ʼ״̬����KF-EKF����õ�Ŀ������״̬
%   ���: 
%       m_Alpha��״̬�м���̣�[x,y,z,M11,M22,M33,M12,M13,M23] ����ο���֪���������5.2.2
%   ����: 
%       Scene���������������в����Լ����շ�������
%           Scene.result.BxyzΪ���������ݣ�Scene.result.Bxyz0Ϊԭʼ���ݡ�
%       Alpha0��Ŀ���ʼ״̬

% 
% Mt=Scene.model.metal.M;
% AlphaError=[0;0;0  ;0;0;0;0;0;0];
% Alpha0=[Scene.model.metal.postion;Mt(1,1);Mt(2,2);Mt(3,3);Mt(1,2);Mt(1,3);Mt(2,3)]+AlphaError;
IterCount=length(Scene.dataconf.v_y)*length(Scene.dataconf.v_x);
m_Alpha=zeros(IterCount+1,length(Alpha0));
k=1;
m_Alpha(k,:)=Alpha0';

R0=1;        % �۲����ֲ�����
Q0=10^-8;            % Ԥ�����ֲ�����
P0=10;           % ��ʼ����
Line_Pk=P0*eye(6);      % ���Բ���Э�������
NonLine_Pk=P0*eye(9);   % �����Բ���Э�������
%% EKF
Peak_amp=max(max(abs(Scene.result.nHz)));
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<Peak_amp/8.2
                  continue;
              end
             
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %��ǰ��Ȧλ��
              curAlpha=m_Alpha(k,:)';
               k=k+1;
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %ʵ�ʹ۲�����
            %% EKF���Ѵż���������v_M������֪������ģ���з����Բ�������v_r
              xk=curAlpha;                                           % Ԥ��λʸr
              NonLine_Pk=Q0*eye(length(xk))+NonLine_Pk;         % Э�������Ԥ��
              [preZ,Hk]=jaccsd_M(Scene.model.detector,Pk,xk,Scene);   %lw �޸� ����²��� Scene����������
%               Z=[real(Z);imag(Z)];
%               Hk=[real(Hk);imag(Hk)];
%               preZ=[real(preZ);imag(preZ)];
              Z=[real(Z);];
              Hk=[real(Hk);];
              preZ=[real(preZ);];
              K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+R0*eye(length(preZ)));    % ����������
              xk=xk+K*(Z-preZ);                                  % ����xk
              NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;         % ����Э�������
%               M=eye(length(NonLine_Pk))-K*Hk;
%               NonLine_Pk=M*NonLine_Pk*M'+K*R0*eye(length(Z))*K'; % ����Э�������,Լɪ���ȶ���
            %% ��¼���ε�������
             m_Alpha(k,:)=xk;
    end
 end
 m_Alpha=m_Alpha(1:k,:);
% mean(m_Alpha)
% m_Alpha(k,:)


% v_M=GetM(Scene,m_Alpha(k,1:3)');
% M=[v_M(1) v_M(4) v_M(5);...
%    v_M(4) v_M(2) v_M(6);...
%    v_M(5) v_M(6) v_M(3)];
% [V,D]=eig(M)
end

