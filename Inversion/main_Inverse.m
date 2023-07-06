%% 
clc;
clear;
close all;
addpath('../0电磁模型/matlab');     % add model libary
%% 参数定义
% 主场参数
Coil=struct();
Coil.I=20;       % 线圈电流
Coil.R=0.5;     % 线圈半径/m
Coil.f=1000;     % 信号频率
Coil.Postion=[0 0 0];   %线圈位置
% 目标物体(等效的偶极子)参数
Target=struct();
Target.Postion=[0 0 -1];   % 目标位置(x,y,z)（线圈中心为坐标原点）
Target.MagPolar=[ 1 2 3];  % 目标三轴磁极化率,依次x,y,z
Target.Theta=0;            % 俯仰角        
Target.Phi=0;              % 滚转角 ps:所有角度为0时，代表三轴极化方向与线圈坐标系的x,y,z重合
Target.Psi=0;              % 航向角（垂直放置的物体是没有的）

% 观测平面参数
Plane=struct();
Plane.h=0;                % 观测平面高度，即z坐标
Plane.Xrange=[-2 2];      % 观测平面x方向范围
Plane.Xinterval=0.5;      % x方向样点间隔         %lw 间隔滤波里是0.1
Plane.Yrange=[-2 2];      % 观测平面y方向范围
Plane.Yinterval=0.5;      % y方向样点间隔

Plane.v_x=Plane.Xrange(1):Plane.Xinterval:Plane.Xrange(2);  % x方向取点坐标
Plane.v_y=Plane.Yrange(1):Plane.Yinterval:Plane.Yrange(2);  % y方向取点坐标
[Plane.m_x,Plane.m_y]=meshgrid(Plane.v_x,Plane.v_y);        % 观测面扫描网格

%反演参数
Hx=zeros(size(Plane.m_x));
Hy=zeros(size(Plane.m_x));
Hz=zeros(size(Plane.m_x));
k=1;                % 数组组数索引
Coils(length(Plane.v_y)*length(Plane.v_x))=struct(Coil);

%% 2/D区域扫描收集二次场数据
    for i=1:length(Plane.v_y)          % 对二维区域的单元格进行扫描，相对于移动线圈
        for j=1:length(Plane.v_x)
            x=Plane.v_x(j);            % 当前线圈的x坐标
            y=Plane.v_y(i);            % 当前线圈的y坐标
            Coil.Postion(1:2)=[x y ];  % 线圈位置，默认高度为0保持不变   
            
            [Hx(i,j),Hy(i,j),Hz(i,j)]=SecondField(Coil,Target,struct('Postion', Coil.Postion));           % 计算当前位置二次场
            if abs(Hz(i,j))>1
                Coils(k)=Coil;
                H(k)=abs(Hz(i,j));
                k=k+1;
            end
            
        end
    end
    Coils=Coils(1:k-1);
    
   Ht=sqrt(abs(Hx).^2+abs(Hy).^2+abs(Hz).^2);
  
mu=0;
sigma=1;        %lw 滤波是10
nHx=Hx+normrnd(mu,sigma,size(Hx));
nHy=Hy+normrnd(mu,sigma,size(Hx));
nHz=Hz+normrnd(mu,sigma,size(Hx));
nHt=sqrt(nHz.^2+nHx.^2+nHy.^2);
    %% 成图结果观测
figure('name','3-D figure of orignal data');
subplot(2,2,1);
surf(Plane.m_x,Plane.m_y,abs(nHt(:,:,1)));
title('orignal data');
xlabel('x/m');
ylabel('y/m');
zlabel('B/nT');
subplot(2,2,2);
surf(Plane.m_x,Plane.m_y,abs(nHx(:,:,1)));
title('orignal data');
xlabel('x/m');
ylabel('y/m');
zlabel('B/nT');
subplot(2,2,3);
surf(Plane.m_x,Plane.m_y,abs(nHy(:,:,1)));
title('orignal data');
xlabel('x/m');
ylabel('y/m');
zlabel('B/nT');
subplot(2,2,4);
surf(Plane.m_x,Plane.m_y,abs(nHz(:,:,1)));
title('orignal data');
xlabel('x/m');
ylabel('y/m');
zlabel('B/nT');
%% 反演
    %% 迭代参数
     IterTimes=20;         % 迭代次数 
     Mu=4*10^13;                % 步长因子
     AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % 真值
     Alpha0=[0.2 0.4 -1.5 0 0 0 1 2 3];   % 迭代初值，依次为x,y,z,θ，φ，ψ，βx,βy,βz
     m_Alpha=zeros(IterTimes+1,length(Alpha0));      % 存储迭代过程
     InvHz=zeros(size(Hz));                         % 反演过程中生成的扫描场(当前只考虑z方向)
     HzGrad=zeros(IterTimes,length(Alpha0));        % 代价函数的梯度
     v_CostValue=zeros(IterTimes,1);                % 代价函数在迭代过程中的值
     m_Alpha(1,:)=Alpha0;
    

     %% 梯度下降法
     
%      for k=1:IterTimes
%         for i=1:length(Plane.v_y)          % 根据当前迭代参数Alpha计算得到场的分布，以及梯度
%             for j=1:length(Plane.v_x)
%                  x=Plane.v_x(j);            % 当前扫描平面的x坐标
%                  y=Plane.v_y(i);            % 当前扫描平面的y坐标
%                  Coil.Postion(1:2)=[x y ];  % 调整线圈位置，默认高度为0保持不变  
%                  [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,m_Alpha(k,:)); % 求H当前位置场的分布
%                  [temp1,temp2,Grad]=HGradAlpha(Coil,m_Alpha(k,:));   % 求当前位置的梯度
%                  e=abs(InvHz(i,j))-abs(Hz(i,j));
%                  HzGrad(k,:)=HzGrad(k,:)+2*e*abs(Grad);
%                  v_CostValue(k)=v_CostValue(k)+e^2;
%                  
%             end
%         end
%         HzGrad(k,:)=HzGrad(k,:)/(i*j);
%         m_Alpha(k+1,:)=m_Alpha(k,:)-Mu.*HzGrad(k,:);
%         
%      end
     %% 牛顿迭代法（BFGS）
%      [Alpha,m_Alpha,HzGrad,v_CostValue]=BFGS(Alpha0,Hz,Coil,Plane,0.01,10^-30,IterTimes);
%      [mincost index]=min(v_CostValue);
%      Alpha=m_Alpha(index,:);

     %% 牛顿迭代法(Levenberg-Marquardt)
     H=awgn(H,10,'measured');
	fun = @(x)HzAlpha(x,Coils)-H';
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
[x,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(fun,Alpha0,[],[],options);
resnorm
output
x
AlphaTrue
error=abs(AlphaTrue-x)
