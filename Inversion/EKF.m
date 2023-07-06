function [m_Alpha] = EKF(Scene,Alpha0)
%给定仿真数据和探测目标初始状态，用KF-EKF计算得到目标最终状态
%   输出: 
%       m_Alpha，状态中间过程，[x,y,z,M11,M22,M33,M12,M13,M23] 具体参考刘知洋毕设论文5.2.2
%   输入: 
%       Scene，包括仿真中所有参数以及最终仿真数据
%           Scene.result.Bxyz为带噪声数据，Scene.result.Bxyz0为原始数据。
%       Alpha0，目标初始状态

% 
% Mt=Scene.model.metal.M;
% AlphaError=[0;0;0  ;0;0;0;0;0;0];
% Alpha0=[Scene.model.metal.postion;Mt(1,1);Mt(2,2);Mt(3,3);Mt(1,2);Mt(1,3);Mt(2,3)]+AlphaError;
IterCount=length(Scene.dataconf.v_y)*length(Scene.dataconf.v_x);
m_Alpha=zeros(IterCount+1,length(Alpha0));
k=1;
m_Alpha(k,:)=Alpha0';

R0=1;        % 观测误差分布方差
Q0=10^-8;            % 预测误差分布方差
P0=10;           % 初始方差
Line_Pk=P0*eye(6);      % 线性部分协方差矩阵
NonLine_Pk=P0*eye(9);   % 非线性部分协方差矩阵
%% EKF
Peak_amp=max(max(abs(Scene.result.nHz)));
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<Peak_amp/8.2
                  continue;
              end
             
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %当前线圈位置
              curAlpha=m_Alpha(k,:)';
               k=k+1;
              Z=[Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];   %实际观测向量
            %% EKF：把磁极化率向量v_M当做已知，估计模型中非线性参数部分v_r
              xk=curAlpha;                                           % 预测位矢r
              NonLine_Pk=Q0*eye(length(xk))+NonLine_Pk;         % 协方差矩阵预测
              [preZ,Hk]=jaccsd_M(Scene.model.detector,Pk,xk,Scene);   %lw 修改 添加新参数 Scene！！！！！
%               Z=[real(Z);imag(Z)];
%               Hk=[real(Hk);imag(Hk)];
%               preZ=[real(preZ);imag(preZ)];
              Z=[real(Z);];
              Hk=[real(Hk);];
              preZ=[real(preZ);];
              K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+R0*eye(length(preZ)));    % 卡尔曼增益
              xk=xk+K*(Z-preZ);                                  % 更新xk
              NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;         % 更新协方差矩阵
%               M=eye(length(NonLine_Pk))-K*Hk;
%               NonLine_Pk=M*NonLine_Pk*M'+K*R0*eye(length(Z))*K'; % 更新协方差矩阵,约瑟夫稳定化
            %% 记录本次迭代参数
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

