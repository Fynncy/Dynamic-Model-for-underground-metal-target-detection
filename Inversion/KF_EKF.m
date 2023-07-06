function [m_Alpha] = KF_EKF(Scene,Alpha0)
%给定仿真数据和探测目标初始状态，用KF-EKF计算得到目标最终状态
%   输出: 
%       m_Alpha，状态中间过程，[x,y,z,M11,M22,M33,M12,M13,M23] 具体参考刘知洋毕设论文5.2.2
%   输入: 
%       Scene，包括仿真中所有参数以及最终仿真数据
%           Scene.result.Bxyz为带噪声数据，Scene.result.Bxyz0为原始数据。
%       Alpha0，目标初始状态


% Mt=Scene.model.metal.M;
% AlphaError=[0;0;0  ;0;0.02;0.03;0;0;0];
% Alpha0=[Scene.model.metal.postion;Mt(1,1);Mt(2,2);Mt(3,3);Mt(1,2);Mt(1,3);Mt(2,3)]+AlphaError;
IterCount=length(Scene.dataconf.v_y)*length(Scene.dataconf.v_x);
m_Alpha=zeros(IterCount+1,length(Alpha0));
k=1;
m_Alpha(k,:)=Alpha0';

R0=10;        % 观测误差分布方差
Q0=10^-8;            % 预测误差分布方差
%Q0=0;
P0=1;           % 初始方差
Line_Pk=P0*eye(6);      % 线性部分协方差矩阵
NonLine_Pk=P0*eye(3);   % 非线性部分协方差矩阵
%% EK-EKF
Peak_amp=max(max(abs(Scene.result.nHz)));
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<Peak_amp/8.2
                  continue;
              end
             
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %当前线圈位置
              curAlpha=m_Alpha(k,:)';
               k=k+1;
            %% EKF：把磁极化率向量v_M当做已知，估计模型中非线性参数部分v_r
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %实际观测向量
              v_r=curAlpha(1:3);
              v_M=curAlpha(4:9);
              xk=v_r;                                           % 预测位矢r
              NonLine_Pk=Q0*eye(length(xk))+NonLine_Pk;         % 协方差矩阵预测
              [preZ,Hk]=jaccsd_M(Scene.model.detector,Pk,curAlpha,Scene.model.detector.theta(i,j),Scene.model.detector.phi(i,j));
              Hk=Hk(:,1:3);                                     % 仅取关于位置的雅格比行列式
              Z=real(Z);
              Hk=real(Hk);
              preZ=real(preZ);
              K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+R0*eye(length(preZ)));    % 卡尔曼增益
              xk=xk+K*(Z-preZ);                                  % 更新xk
              NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;             % 更新协方差矩阵
            %% KF：把当前位矢r作为已知，估计模型中线性参数部分v_M 
              v_r=xk;
              [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3),Scene.model.detector.theta(i,j),Scene.model.detector.phi(i,j));      %计算格式函数Gk
              %Bp=FirstField(Scene.model.detector.I,Scene.model.detector.R,Scene.model.detector.F,v_r-Pk);   % 计算主场Bp
              Bp=FirstField_wz_matrix(Scene.model.detector.I,Scene.model.detector.R,(v_r-Pk)',Scene.model.detector.theta(i,j),Scene.model.detector.phi(i,j));
              Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...   % 计算Wk
                  0 Bp(2) 0 Bp(1) 0 Bp(3);...
                  0 0 Bp(3) 0 Bp(1) Bp(2)];
              Hk=Gk*Wk;                      % 观测模型中的转移矩阵
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %实际观测向量
              xk=v_M;                                      %预测
              Line_Pk=Q0*eye(length(Line_Pk))+Line_Pk;     % pk|k-1=Qk-1+Fk*Pk-1|k-1*Fk'，预测协方差矩阵
              K=Line_Pk*Hk'/(Hk*Line_Pk*Hk'+R0*eye(3));    % 卡尔曼增益
              xk=xk+K*(Z-Hk*xk);                           % 更新xk
              Line_Pk=Line_Pk-K*Hk*Line_Pk;                % 更新方差矩阵
              v_M=xk;
            %% 记录本次迭代参数
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

