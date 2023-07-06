function [m_Alpha] = KF_EKF_lw(Scene,Alpha0)
%�����������ݺ�̽��Ŀ���ʼ״̬����KF-EKF����õ�Ŀ������״̬
%   ���: 
%       m_Alpha��״̬�м���̣�[x,y,z,M11,M22,M33,M12,M13,M23] ����ο���֪���������5.2.2
%   ����: 
%       Scene���������������в����Լ����շ�������
%           Scene.result.BxyzΪ���������ݣ�Scene.result.Bxyz0Ϊԭʼ���ݡ�
%       Alpha0��Ŀ���ʼ״̬


% Mt=Scene.model.metal.M;
% AlphaError=[0;0;0  ;0;0.02;0.03;0;0;0];
% Alpha0=[Scene.model.metal.postion;Mt(1,1);Mt(2,2);Mt(3,3);Mt(1,2);Mt(1,3);Mt(2,3)]+AlphaError;
IterCount=length(Scene.dataconf.v_y)*length(Scene.dataconf.v_x);
m_Alpha=zeros(IterCount+1,length(Alpha0));
k=1;
m_Alpha(k,:)=Alpha0';

%R0=1;
Q0=10^-8;

% Rv_r=10^-4;        % �۲����ֲ�����
Rv_r=[1,0,0; 0,1,0; 0,0,1];
% Qv_r=[10^-8,0,0; 0,10^-8,0; 0,0,10^-8];            % Ԥ�����ֲ�����

Rv_M=[1,0,0; 0,1,0; 0,0,1];
% Qv_M=[10^-8,0,0,0,0,0; 0,10^-8,0,0,0,0; 0,0,10^-8,0,0,0; 0,0,0,10^-8,0,0; 0,0,0,0,10^-8,0; 0,0,0,0,0,10^-8];
% Rv_M=10^-5;
% Qv_M=10^-8; 

P0=10;           % ��ʼ����
% P1=10;
Line_Pk=P0*eye(6);      % ���Բ���Э�������     %lw  ���� ���Э�������
NonLine_Pk=P0*eye(3);   % �����Բ���Э�������   %lw ������ ���Э�������
%% EK-EKF
Peak_amp=max(max(abs(Scene.result.nHz)));
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<Peak_amp/8.2
                  continue;
              end
             
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %��ǰ��Ȧλ��
              curAlpha=m_Alpha(k,:)';
              k=k+1;
            %% EKF���Ѵż���������v_M������֪������ģ���з����Բ�������v_r
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %ʵ�ʹ۲�����
              v_r=curAlpha(1:3);
              v_M=curAlpha(4:9);
              xk=v_r;                                           % Ԥ��λʸr
              NonLine_Pk=Q0+NonLine_Pk;         % Э�������Ԥ�� %lw lzy 4-30 ԭQv_r*eye(length(xk))
              [preZ,Hk]=jaccsd_M(Scene.model.detector,Pk,curAlpha,Scene);   %lw  �޸� ����²���Scene %lw lzy 4-34
              Hk=Hk(:,1:3);                                     % ��ȡ����λ�õ��Ÿ������ʽ
              Z=real(Z);
              Hk=real(Hk);
              preZ=real(preZ);
              K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+Rv_r);    % ���������� %lw lzy 4-33 ԭRv_r*eye(length(preZ))
              xk=xk+K*(Z-preZ);                                  % ����xk %lw lzy 4-31
              NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;             % ����Э������� %lw lzy 4-32
            %% KF���ѵ�ǰλʸr��Ϊ��֪������ģ�������Բ�������v_M 
              v_r=xk;     %lw  lzy v_r��֪
              [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3),Scene.model.detector.theta,Scene.model.detector.phi);      %�����ʽ����Gk %lw�޸ģ������̬�ǲ���
              Bp=FirstField(Scene.model.detector.I,Scene.model.detector.R,Scene.model.detector.F,v_r-Pk);   % ��������Bp
              Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...   % ����Wk
                  0 Bp(2) 0 Bp(1) 0 Bp(3);...
                  0 0 Bp(3) 0 Bp(1) Bp(2)];
              Hk=Gk*Wk;                      % �۲�ģ���е�ת�ƾ���  %lw lzy 4-21
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %ʵ�ʹ۲�����
              xk=v_M;                                      %Ԥ��
              Line_Pk=Q0+Line_Pk;     % pk|k-1=Qk-1+Fk*Pk-1|k-1*Fk'��Ԥ��Э�������  %lw lzy 4-23 ԭQv_M*eye(length(Line_Pk))              
              K=Line_Pk*Hk'/(Hk*Line_Pk*Hk'+Rv_M);    % ���������� %lw lzy 4-26 ����Hk=Gk*Wk ԭRv_M*eye(3)
              xk=xk+K*(Z-Hk*xk);                           % ����xk
              Line_Pk=Line_Pk-K*Hk*Line_Pk;                % ���·������
              v_M=xk;
            %% ��¼���ε�������
             m_Alpha(k,:)=[v_r' v_M'];
    end
 end
 m_Alpha=m_Alpha(1:k,:);
%  k
% mean(m_Alpha)
% m_Alpha(k,:)
% v_M=m_Alpha(k,4:9);
% M=[v_M(1) v_M(4) v_M(5);...
%    v_M(4) v_M(2) v_M(6);...
%    v_M(5) v_M(6) v_M(3)];
% [V,D]=eig(M)
end

