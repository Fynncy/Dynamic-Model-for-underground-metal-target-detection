%% 
clc;
clear;
close all;
addpath('../0���ģ��/matlab');     % add model libary
%% ��������
% ��������
Coil=struct();
Coil.I=20;       % ��Ȧ����
Coil.R=0.5;     % ��Ȧ�뾶/m
Coil.f=1000;     % �ź�Ƶ��
Coil.Postion=[0 0 0];   %��Ȧλ��
% Ŀ������(��Ч��ż����)����
Target=struct();
Target.Postion=[0 0 -1];   % Ŀ��λ��(x,y,z)����Ȧ����Ϊ����ԭ�㣩
Target.MagPolar=[ 1 2 3];  % Ŀ������ż�����,����x,y,z
Target.Theta=0;            % ������        
Target.Phi=0;              % ��ת�� ps:���нǶ�Ϊ0ʱ���������Ἣ����������Ȧ����ϵ��x,y,z�غ�
Target.Psi=0;              % ����ǣ���ֱ���õ�������û�еģ�

% �۲�ƽ�����
Plane=struct();
Plane.h=0;                % �۲�ƽ��߶ȣ���z����
Plane.Xrange=[-2 2];      % �۲�ƽ��x����Χ
Plane.Xinterval=0.5;      % x����������         %lw ����˲�����0.1
Plane.Yrange=[-2 2];      % �۲�ƽ��y����Χ
Plane.Yinterval=0.5;      % y����������

Plane.v_x=Plane.Xrange(1):Plane.Xinterval:Plane.Xrange(2);  % x����ȡ������
Plane.v_y=Plane.Yrange(1):Plane.Yinterval:Plane.Yrange(2);  % y����ȡ������
[Plane.m_x,Plane.m_y]=meshgrid(Plane.v_x,Plane.v_y);        % �۲���ɨ������

%���ݲ���
Hx=zeros(size(Plane.m_x));
Hy=zeros(size(Plane.m_x));
Hz=zeros(size(Plane.m_x));
k=1;                % ������������
Coils(length(Plane.v_y)*length(Plane.v_x))=struct(Coil);

%% 2/D����ɨ���ռ����γ�����
    for i=1:length(Plane.v_y)          % �Զ�ά����ĵ�Ԫ�����ɨ�裬������ƶ���Ȧ
        for j=1:length(Plane.v_x)
            x=Plane.v_x(j);            % ��ǰ��Ȧ��x����
            y=Plane.v_y(i);            % ��ǰ��Ȧ��y����
            Coil.Postion(1:2)=[x y ];  % ��Ȧλ�ã�Ĭ�ϸ߶�Ϊ0���ֲ���   
            
            [Hx(i,j),Hy(i,j),Hz(i,j)]=SecondField(Coil,Target,struct('Postion', Coil.Postion));           % ���㵱ǰλ�ö��γ�
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
sigma=1;        %lw �˲���10
nHx=Hx+normrnd(mu,sigma,size(Hx));
nHy=Hy+normrnd(mu,sigma,size(Hx));
nHz=Hz+normrnd(mu,sigma,size(Hx));
nHt=sqrt(nHz.^2+nHx.^2+nHy.^2);
    %% ��ͼ����۲�
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
%% ����
    %% ��������
     IterTimes=20;         % �������� 
     Mu=4*10^13;                % ��������
     AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % ��ֵ
     Alpha0=[0.2 0.4 -1.5 0 0 0 1 2 3];   % ������ֵ������Ϊx,y,z,�ȣ��գ��ף���x,��y,��z
     m_Alpha=zeros(IterTimes+1,length(Alpha0));      % �洢��������
     InvHz=zeros(size(Hz));                         % ���ݹ��������ɵ�ɨ�賡(��ǰֻ����z����)
     HzGrad=zeros(IterTimes,length(Alpha0));        % ���ۺ������ݶ�
     v_CostValue=zeros(IterTimes,1);                % ���ۺ����ڵ��������е�ֵ
     m_Alpha(1,:)=Alpha0;
    

     %% �ݶ��½���
     
%      for k=1:IterTimes
%         for i=1:length(Plane.v_y)          % ���ݵ�ǰ��������Alpha����õ����ķֲ����Լ��ݶ�
%             for j=1:length(Plane.v_x)
%                  x=Plane.v_x(j);            % ��ǰɨ��ƽ���x����
%                  y=Plane.v_y(i);            % ��ǰɨ��ƽ���y����
%                  Coil.Postion(1:2)=[x y ];  % ������Ȧλ�ã�Ĭ�ϸ߶�Ϊ0���ֲ���  
%                  [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,m_Alpha(k,:)); % ��H��ǰλ�ó��ķֲ�
%                  [temp1,temp2,Grad]=HGradAlpha(Coil,m_Alpha(k,:));   % ��ǰλ�õ��ݶ�
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
     %% ţ�ٵ�������BFGS��
%      [Alpha,m_Alpha,HzGrad,v_CostValue]=BFGS(Alpha0,Hz,Coil,Plane,0.01,10^-30,IterTimes);
%      [mincost index]=min(v_CostValue);
%      Alpha=m_Alpha(index,:);

     %% ţ�ٵ�����(Levenberg-Marquardt)
     H=awgn(H,10,'measured');
	fun = @(x)HzAlpha(x,Coils)-H';
    options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
[x,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(fun,Alpha0,[],[],options);
resnorm
output
x
AlphaTrue
error=abs(AlphaTrue-x)
