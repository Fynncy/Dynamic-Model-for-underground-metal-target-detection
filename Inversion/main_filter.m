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
Plane.Xinterval=0.1;      % x方向样点间隔
Plane.Yrange=[-2 2];      % 观测平面y方向范围
Plane.Yinterval=0.1;      % y方向样点间隔

Plane.v_x=Plane.Xrange(1):Plane.Xinterval:Plane.Xrange(2);  % x方向取点坐标
Plane.v_y=Plane.Yrange(1):Plane.Yinterval:Plane.Yrange(2);  % y方向取点坐标
[Plane.m_x,Plane.m_y]=meshgrid(Plane.v_x,Plane.v_y);        % 观测面扫描网格

%反演参数


 
Hx=zeros(size(Plane.m_x));
Hy=zeros(size(Plane.m_x));
Hz=zeros(size(Plane.m_x));
k=0;                % 数组组数索引
Coils(length(Plane.v_y)*length(Plane.v_x))=struct(Coil);
Target.Dipole=CalMoment(Target.MagPolar,FirstField(Coil.I,Coil.R,Coil.f,Target.Postion-Coil.Postion),Target.Theta,Target.Phi,Target.Psi);
%% 2/D区域扫描收集二次场数据
    for i=1:length(Plane.v_y)          % 对二维区域的单元格进行扫描，相对于移动线圈
        for j=1:length(Plane.v_x)
            k=k+1;
            x=Plane.v_x(j);            % 当前线圈的x坐标
            y=Plane.v_y(i);            % 当前线圈的y坐标
            Coil.Postion(1:2)=[x y ];  % 改变发现线圈位置，使其和扫描位置一致。否则主场不发生变化，偶极子模型。  
            RcPostion=[x y Plane.h];    %接收线圈位置
            Coils(k)=Coil;
            [Hx(i,j),Hy(i,j),Hz(i,j)]=SecondField(Coil,Target,struct('Postion', RcPostion));           % 计算当前位置二次场
          % H(k)=abs(Hz(i,j));
        end
    end
   Ht=sqrt(abs(Hx).^2+abs(Hy).^2+abs(Hz).^2);
%% 噪声
mu=0;
sigma=10;    % 噪声分布方差
clear j;
nHx=Hx+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);
nHy=Hy+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);
nHz=Hz+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);

% nHx=Hx+normrnd(mu,sigma,size(Hx));
% nHy=Hy+normrnd(mu,sigma,size(Hx));
% nHz=Hz+normrnd(mu,sigma,size(Hx));

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
% figure('name','Image');
% title('orignal image of amplitude');
% image(Plane.v_x,Plane.v_y,uint8(mapminmax(Ht(:,:,1),0,255)));
% xlabel('x/m');
% ylabel('y/m');


%% 滤波估值   

    %% KF

    %% EXF
    
    %% KF-EKF
    IterCount=length(Plane.v_y)*length(Plane.v_x);
    AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % 真值
    Alpha0=[0.2 -0.1 -1 0 0 0 1 2 2];   % 迭代初值，依次为x,y,z,θ，φ，ψ，βx,βy,βz
    R0=100;        % 观测误差分布方差
    Q0=1;            % 预测误差分布方差
    P0=1;           % 初始方差
    m_Alpha=zeros(IterCount+1,length(Alpha0));  
    k=1;
    m_Alpha(k,:)=Alpha0;
    Line_Pk=P0*eye(3);      % 线性部分协方差矩阵
    NonLine_Pk=P0*eye(6);   % 非线性部分协方差矩阵
     for i=1:length(Plane.v_y)          
        for j=1:length(Plane.v_x)       % 遍历每一个发射线圈的位置
            if abs(nHt(i,j))<100
                continue;
            end
           curX=Plane.v_x(j);
           curY=Plane.v_y(i);
           TxPostion=[ curX curY Coil.Postion(3)]'; % 发射线圈位置
           RxPostion=[ curX curY Plane.h]';         % 接收线圈位置
           curAlpha=m_Alpha(k,:);                   % 当前状态向量
           k=k+1;
           %% KF：通过估计磁矩[mx,my,mz]得到极化率[βx,βy,βz]
           v_targetPostion=curAlpha(1:3);           % 取当前物体位置矢量
           v_targetPostrue=curAlpha(4:6);           % 取当前物体姿态
           v_targetMagPolar=curAlpha(7:9);          % 去当前物体磁极化率
           v_Br=FirstField(Coil.I,Coil.R,Coil.f,v_targetPostion'-TxPostion);   % 计算在物体位置处的主场
           % 计算感应磁偶极矩
           v_inducedDipMoment=CalMoment(v_targetMagPolar,v_Br,v_targetPostrue(1),v_targetPostrue(2),v_targetPostrue(3)); 
           % 预测-更新磁偶极矩
           xk=v_inducedDipMoment;                       % 预测xk|k-1
           [bx,by,bz,Hk]=HFieldModel([v_targetPostion(:) v_inducedDipMoment(:)],RxPostion(1),RxPostion(2),RxPostion(3),v_targetPostrue(1),v_targetPostrue(2));    % 计算hk以及Hk（即根据磁矩和位置矢量计算观测场）
           Line_Pk=Q0*eye(length(Line_Pk))+Line_Pk;     % pk|k-1=Qk-1+Fk*Pk-1|k-1*Fk'，预测协方差矩阵 %lw 4-23
           K=Line_Pk*Hk'/(Hk*Line_Pk*Hk'+R0*eye(3));    % 卡尔曼增益  %lw 4-26
           Z=[nHx(i,j);nHy(i,j);nHz(i,j)];
           xk=xk+K*(Z-[bx;by;bz]);        %更新xk   %lw 4-24
           Line_Pk=Line_Pk-K*Hk*Line_Pk;                % 更新方差矩阵 %lw 4-25
           v_inducedDipMoment=xk;
           Rt=RotationTensor(v_targetPostrue(1),v_targetPostrue(2),v_targetPostrue(3));
           temp1=Rt'*v_inducedDipMoment;
           temp2=Rt'*v_Br';
           for p=1:length(temp2) % 计算极化率βx,βy,βz
               if abs(temp2(p))>2
                   v_targetMagPolar(p)=temp1(p)/temp2(p);
               end
           end
           %% EKF：估计[x,y,z,θ，φ，ψ]
           xk=[v_targetPostion v_targetPostrue]';           % 预测
           NonLine_Pk=Q0*eye(length(xk))+NonLine_Pk;        % 协方差矩阵预测 %lw 4-30
           Coil.Postion=TxPostion;
           [Z,Hk]=jaccsd(Coil,xk',v_targetMagPolar);
           hk=[real(nHx(i,j));real(nHy(i,j));real(nHz(i,j))];%imag(nHx(i,j));imag(nHy(i,j));imag(nHz(i,j))];
           Z=[real(Z)];
           Hk=[real(Hk)];
           K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+R0*eye(length(Z)));    % 卡尔曼增益  %lw 4-33
           xk=xk+K*(hk-Z);                                  % 更新xk
           NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;         % 更新协方差矩阵  %lw 4-32
%            M=eye(length(NonLine_Pk))-K*Hk;
%            NonLine_Pk=M*NonLine_Pk*M'+K*R0*eye(length(Z))*K'; % 更新协方差矩阵,约瑟夫稳定化
           
           v_targetPostion=xk(1:3)';
           v_targetPostrue=xk(4:6)';
          %% 记录本次迭代参数
          
          m_Alpha(k,:)=[v_targetPostion v_targetPostrue v_targetMagPolar];
        end
     end
    m_Alpha=m_Alpha(1:k,:);
    mean(m_Alpha)
    m_Alpha(end,:)
%% 反演
    %% 迭代参数
%      IterTimes=20;         % 迭代次数 
%      Mu=4*10^13;                % 步长因子
%      AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % 真值
%      Alpha0=[0.3 0.2 -1 0 0 0 0.8 1 1];   % 迭代初值，依次为x,y,z,θ，φ，ψ，βx,βy,βz
%      m_Alpha=zeros(IterTimes+1,length(Alpha0));      % 存储迭代过程
%      InvHz=zeros(size(Hz));                         % 反演过程中生成的扫描场(当前只考虑z方向)
%      HzGrad=zeros(IterTimes,length(Alpha0));        % 代价函数的梯度
%      v_CostValue=zeros(IterTimes,1);                % 代价函数在迭代过程中的值
%      m_Alpha(1,:)=Alpha0;
    

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
%      H=awgn(H,2000,'measured');
% 	fun = @(x)HzAlpha(x,Coils)-H';
%     options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
% [x,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(fun,Alpha0,[],[],options);
% resnorm
% output
% x
% AlphaTrue
% error=abs(AlphaTrue-x)
