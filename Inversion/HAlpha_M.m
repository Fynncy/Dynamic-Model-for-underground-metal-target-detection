function [B] = HAlpha_M(detector,Pk,Alpha,theta,phi)
%给定线圈参数，根据模型参数计算二次场 （具体理论参考刘知洋毕设论文-2.3节）
%   输出: 
%       B，二次场
%   输入: 
%       detector，发射器参数
%       pk,发射器位置
%       Alpha，1*9向量，依次为x,y,z，M11,M22,M33,M12,M13,M23 

  v_r=Alpha(1:3);
  v_M=Alpha(4:9);
  v_r=v_r(:);
  v_M=v_M(:);
  Pk=Pk(:);
  [~,~,~,Gk]=HFieldModel([v_r v_r],Pk(1),Pk(2),Pk(3),theta,phi);      %计算格式函数Gk
  %Bp=FirstField(detector.I,detector.R,detector.F,v_r-Pk);   % 计算主场Bp %lw
  %注释掉了变成下行：
  Bp=FirstField_wz_matrix(detector.I,detector.R,(v_r-Pk)',theta,phi);
  Wk=[Bp(1) 0 0 Bp(2) Bp(3) 0;...   % 计算Wk
      0 Bp(2) 0 Bp(1) 0 Bp(3);...
      0 0 Bp(3) 0 Bp(1) Bp(2)];
   B=Gk*Wk*v_M;

end

