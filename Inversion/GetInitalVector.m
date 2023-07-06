function [v_inital] = GetInitalVector(Scene)
%给定仿真数据，计算探测目标初始状态
%   输出: 
%       v_inital，初始状态，[x,y,z,M11,M22,M33,M12,M13,M23] 具体参考论文
%   输入: 
%       Scene，包括仿真中所有参数以及最终仿真数据
%       Scene.result.Bxyz为带噪声数据，Scene.result.Bxyz0为原始数据。
%       
v_inital=zeros(9,1);
%% 估计水平位置
x=sum(sum(Scene.result.Bxyz.*Scene.dataconf.m_x))/sum(sum(Scene.result.Bxyz));
y=sum(sum(Scene.result.Bxyz.*Scene.dataconf.m_y))/sum(sum(Scene.result.Bxyz));

v_inital(1:2)=[x y]';
%% 估计深度
area=Scene.dataconf.linespace*Scene.dataconf.interval;
z=-2*sqrt(FWHM(Scene.result.Bxyz,area)/pi);
v_inital(3)=z;
%% 估计磁极化率张量M
Scene.model.metal.M;
r=[x y z]'; %%初始位置矢量由空间响应一阶矩和半峰波宽确定
% r=Scene.model.metal.postion(:); % 初始值为目标真实位置
Bs=sqrt(Scene.result.nHx.^2+Scene.result.nHy.^2+Scene.result.nHz.^2);
W=[];
Y=[];
k=1;
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<100
                  continue;
              end
              k=k+1;
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %当前线圈位置
              [~,~,~,Gk]=HFieldModel([r r],Pk(1),Pk(2),Pk(3));
              Bp=FirstField(Scene.model.detector.I,Scene.model.detector.R,Scene.model.detector.F,r-Pk);
              Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...
                  0 Bp(2) 0 Bp(1) 0 Bp(3);...
                  0 0 Bp(3) 0 Bp(1) Bp(2)];
              W=[W;Gk*Wk];
              Y=[Y;Scene.result.nHx(i,j);Scene.result.nHy(i,j);Scene.result.nHz(i,j);];
              
    end
 end
 v_M=(W'*W)\(W'*Y);
%  k
%  cond(W'*W)
%  det(W'*W)
%  size(W)
%  cond(W'*W)
%  M=[v_M(1) v_M(4) v_M(5);...
%      v_M(4) v_M(2) v_M(6);...
%      v_M(5) v_M(6) v_M(3)];
%  M
%  Scene.model.metal.M
%  M- Scene.model.metal.M
 v_inital(4:9)=v_M(1:6);
end

