function [v_M] = GetM(Scene,r)
%给定仿真数据，计算探测目标初始状态
%   输出: 
%       v_M，磁极化率张量独立元素张量，[M11,M22,M33,M12,M13,M23] 具体参考论文
%   输入: 
%       Scene，包括仿真中所有参数以及最终仿真数据
%       Scene.result.Bxyz为带噪声数据，Scene.result.Bxyz0为原始数据。
%       r，当前位置矢量
%       
%% 估计磁极化率张量M
Scene.model.metal.M;
% r=Scene.model.metal.postion(:); % 初始值为目标真实位置
Bs=sqrt(Scene.result.nHx.^2+Scene.result.nHy.^2+Scene.result.nHz.^2);
W=[];
Y=[];
k=1;
 for i=1:length(Scene.dataconf.v_y)         
     for j=1:length(Scene.dataconf.v_x)
              if abs(Scene.result.nHz(i,j))<10
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
end

