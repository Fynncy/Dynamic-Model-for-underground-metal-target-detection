function [X,Xs,HzGrad,v_CostValue] = BFGS(Alpha,H,Coil,Plane,mu,precision,N)
%BFGS BFGS法求最小二乘解
%   给定给定步长和指定精度，用BFGS算法求解最小二乘解，
%       输出：
%           X，迭代初值
%           Xs,迭代过程
%       输入：    
%           Alpha，迭代初值
%           Coil,线圈
%           H，观测数据(默认取Hz轴)
%           Plane，观测平面
%           mu,步长因子（暂且用1，没求最优化步长因子）
%           precision，精度
%           N,最大迭代次数


%% 初始化参数
X=Alpha(:);                           % x0，迭代初值
Xs=zeros(N,length(Alpha));            % N*m矩阵，记录迭代过程
D0=eye(length(Alpha));                % D0,迭代近似hessian矩阵的逆
v_CostValue=zeros(N,1);               % 代价函数在迭代过程中的值
HzGrad=zeros(N,length(Alpha));        % 代价函数的梯度
InvHz=zeros(size(H));                 % 反演过程中生成的扫描场(当前只考虑z方向)
%% 计算gk
k=1;
        for i=1:length(Plane.v_y)          
            for j=1:length(Plane.v_x)
                 x=Plane.v_x(j);            % 当前扫描平面的x坐标
                 y=Plane.v_y(i);            % 当前扫描平面的y坐标
                 Coil.Postion(1:2)=[x y ];  % 调整线圈位置，默认高度为0保持不变  
                 [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,Alpha); % 求H当前位置场的分布
                 [temp1,temp2,Grad]=HGradAlpha(Coil,Alpha);   % 求当前位置的梯度
                 e=abs(InvHz(i,j))-abs(H(i,j));
                 HzGrad(k,:)=HzGrad(k,:)+2*e*Grad;
                 v_CostValue(k)=v_CostValue(k)+e^2;
            end
        end
        
gk=HzGrad(k,:)';
deltax=-mu*D0*gk;              % 牛顿方向
dis=sqrt(dot(deltax,deltax));
Xs(1,:)=X';
k=1;
while dis>precision &&k<N  
    k=k+1;
    X=X+deltax;
    Xs(k,:)=X';
        for i=1:length(Plane.v_y)          
            for j=1:length(Plane.v_x)
                 x=Plane.v_x(j);            % 当前扫描平面的x坐标
                 y=Plane.v_y(i);            % 当前扫描平面的y坐标
                 Coil.Postion(1:2)=[x y ];  % 调整线圈位置，默认高度为0保持不变  
                 [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,X'); % 求H当前位置场的分布
                 [temp1,temp2,Grad]=HGradAlpha(Coil,X');   % 求当前位置的梯度
                 e=abs(InvHz(i,j))-abs(H(i,j));
                 HzGrad(k,:)=HzGrad(k,:)+2*e*abs(Grad);
                 v_CostValue(k)=v_CostValue(k)+e^2;
            end
        end
    delta_gk=HzGrad(k,:)'-gk;
    gk=HzGrad(k,:)'; 
    D0=(eye(size(D0))-deltax*delta_gk'/(delta_gk'*deltax))*D0*(eye(size(D0))-delta_gk*deltax'/(delta_gk'*deltax))+deltax*deltax'/(delta_gk'*deltax);
    deltax=-mu*D0*gk;
    dis=sqrt(dot(deltax,deltax));
end
Xs=Xs(1:k,:);
end

