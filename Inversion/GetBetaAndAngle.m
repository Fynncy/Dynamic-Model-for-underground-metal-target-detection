function [Beta,Angles] = GetBetaAndAngle(v_M)
%给定仿真数据，计算探测目标初始状态
%   输出: 
%       v_inital，初始状态，[x,y,z,M11,M22,M33,M12,M13,M23] 具体毕设论文-基于滤波跟踪的电磁探测算法研究
%   输入: 
%       v_M，M11,M22,M33,M12,M13,M23
%       
%% 根据v_M求特征值求解主轴极化率以及姿态角
M=[v_M(1) v_M(4) v_M(5);...
   v_M(4) v_M(2) v_M(6);...
   v_M(5) v_M(6) v_M(3)];
[V,D]=eig(M);                          %lw 求矩阵A的全部特征值，构成对角阵D，并求A的特征向量构成V的列向量
Beta=sort([ D(3,3) D(2,2) D(1,1)]);
Angles(1)=asin(-real(V(1,1)));			%v_M中可能有复数，故特征向量只需取实部
Angles(2)=asin(real(V(2,1))/cos(Angles(1)));
Angles(3)=asin(real(V(1,2))/cos(Angles(2)));
Angles=Angles*180/pi;
Angles(Angles<0)=-1*Angles(Angles<0);		%认为角度只大于0

end

