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
Plane.Xinterval=0.1;      % x����������
Plane.Yrange=[-2 2];      % �۲�ƽ��y����Χ
Plane.Yinterval=0.1;      % y����������

Plane.v_x=Plane.Xrange(1):Plane.Xinterval:Plane.Xrange(2);  % x����ȡ������
Plane.v_y=Plane.Yrange(1):Plane.Yinterval:Plane.Yrange(2);  % y����ȡ������
[Plane.m_x,Plane.m_y]=meshgrid(Plane.v_x,Plane.v_y);        % �۲���ɨ������

%���ݲ���


 
Hx=zeros(size(Plane.m_x));
Hy=zeros(size(Plane.m_x));
Hz=zeros(size(Plane.m_x));
k=0;                % ������������
Coils(length(Plane.v_y)*length(Plane.v_x))=struct(Coil);
Target.Dipole=CalMoment(Target.MagPolar,FirstField(Coil.I,Coil.R,Coil.f,Target.Postion-Coil.Postion),Target.Theta,Target.Phi,Target.Psi);
%% 2/D����ɨ���ռ����γ�����
    for i=1:length(Plane.v_y)          % �Զ�ά����ĵ�Ԫ�����ɨ�裬������ƶ���Ȧ
        for j=1:length(Plane.v_x)
            k=k+1;
            x=Plane.v_x(j);            % ��ǰ��Ȧ��x����
            y=Plane.v_y(i);            % ��ǰ��Ȧ��y����
            Coil.Postion(1:2)=[x y ];  % �ı䷢����Ȧλ�ã�ʹ���ɨ��λ��һ�¡����������������仯��ż����ģ�͡�  
            RcPostion=[x y Plane.h];    %������Ȧλ��
            Coils(k)=Coil;
            [Hx(i,j),Hy(i,j),Hz(i,j)]=SecondField(Coil,Target,struct('Postion', RcPostion));           % ���㵱ǰλ�ö��γ�
          % H(k)=abs(Hz(i,j));
        end
    end
   Ht=sqrt(abs(Hx).^2+abs(Hy).^2+abs(Hz).^2);
%% ����
mu=0;
sigma=10;    % �����ֲ�����
clear j;
nHx=Hx+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);
nHy=Hy+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);
nHz=Hz+normrnd(mu,sigma,size(Hx))*exp(2*pi*rand(1)*j);

% nHx=Hx+normrnd(mu,sigma,size(Hx));
% nHy=Hy+normrnd(mu,sigma,size(Hx));
% nHz=Hz+normrnd(mu,sigma,size(Hx));

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
% figure('name','Image');
% title('orignal image of amplitude');
% image(Plane.v_x,Plane.v_y,uint8(mapminmax(Ht(:,:,1),0,255)));
% xlabel('x/m');
% ylabel('y/m');


%% �˲���ֵ   

    %% KF

    %% EXF
    
    %% KF-EKF
    IterCount=length(Plane.v_y)*length(Plane.v_x);
    AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % ��ֵ
    Alpha0=[0.2 -0.1 -1 0 0 0 1 2 2];   % ������ֵ������Ϊx,y,z,�ȣ��գ��ף���x,��y,��z
    R0=100;        % �۲����ֲ�����
    Q0=1;            % Ԥ�����ֲ�����
    P0=1;           % ��ʼ����
    m_Alpha=zeros(IterCount+1,length(Alpha0));  
    k=1;
    m_Alpha(k,:)=Alpha0;
    Line_Pk=P0*eye(3);      % ���Բ���Э�������
    NonLine_Pk=P0*eye(6);   % �����Բ���Э�������
     for i=1:length(Plane.v_y)          
        for j=1:length(Plane.v_x)       % ����ÿһ��������Ȧ��λ��
            if abs(nHt(i,j))<100
                continue;
            end
           curX=Plane.v_x(j);
           curY=Plane.v_y(i);
           TxPostion=[ curX curY Coil.Postion(3)]'; % ������Ȧλ��
           RxPostion=[ curX curY Plane.h]';         % ������Ȧλ��
           curAlpha=m_Alpha(k,:);                   % ��ǰ״̬����
           k=k+1;
           %% KF��ͨ�����ƴž�[mx,my,mz]�õ�������[��x,��y,��z]
           v_targetPostion=curAlpha(1:3);           % ȡ��ǰ����λ��ʸ��
           v_targetPostrue=curAlpha(4:6);           % ȡ��ǰ������̬
           v_targetMagPolar=curAlpha(7:9);          % ȥ��ǰ����ż�����
           v_Br=FirstField(Coil.I,Coil.R,Coil.f,v_targetPostion'-TxPostion);   % ����������λ�ô�������
           % �����Ӧ��ż����
           v_inducedDipMoment=CalMoment(v_targetMagPolar,v_Br,v_targetPostrue(1),v_targetPostrue(2),v_targetPostrue(3)); 
           % Ԥ��-���´�ż����
           xk=v_inducedDipMoment;                       % Ԥ��xk|k-1
           [bx,by,bz,Hk]=HFieldModel([v_targetPostion(:) v_inducedDipMoment(:)],RxPostion(1),RxPostion(2),RxPostion(3),v_targetPostrue(1),v_targetPostrue(2));    % ����hk�Լ�Hk�������ݴžغ�λ��ʸ������۲ⳡ��
           Line_Pk=Q0*eye(length(Line_Pk))+Line_Pk;     % pk|k-1=Qk-1+Fk*Pk-1|k-1*Fk'��Ԥ��Э������� %lw 4-23
           K=Line_Pk*Hk'/(Hk*Line_Pk*Hk'+R0*eye(3));    % ����������  %lw 4-26
           Z=[nHx(i,j);nHy(i,j);nHz(i,j)];
           xk=xk+K*(Z-[bx;by;bz]);        %����xk   %lw 4-24
           Line_Pk=Line_Pk-K*Hk*Line_Pk;                % ���·������ %lw 4-25
           v_inducedDipMoment=xk;
           Rt=RotationTensor(v_targetPostrue(1),v_targetPostrue(2),v_targetPostrue(3));
           temp1=Rt'*v_inducedDipMoment;
           temp2=Rt'*v_Br';
           for p=1:length(temp2) % ���㼫���ʦ�x,��y,��z
               if abs(temp2(p))>2
                   v_targetMagPolar(p)=temp1(p)/temp2(p);
               end
           end
           %% EKF������[x,y,z,�ȣ��գ���]
           xk=[v_targetPostion v_targetPostrue]';           % Ԥ��
           NonLine_Pk=Q0*eye(length(xk))+NonLine_Pk;        % Э�������Ԥ�� %lw 4-30
           Coil.Postion=TxPostion;
           [Z,Hk]=jaccsd(Coil,xk',v_targetMagPolar);
           hk=[real(nHx(i,j));real(nHy(i,j));real(nHz(i,j))];%imag(nHx(i,j));imag(nHy(i,j));imag(nHz(i,j))];
           Z=[real(Z)];
           Hk=[real(Hk)];
           K=NonLine_Pk*Hk'/(Hk*NonLine_Pk*Hk'+R0*eye(length(Z)));    % ����������  %lw 4-33
           xk=xk+K*(hk-Z);                                  % ����xk
           NonLine_Pk=NonLine_Pk-K*Hk*NonLine_Pk;         % ����Э�������  %lw 4-32
%            M=eye(length(NonLine_Pk))-K*Hk;
%            NonLine_Pk=M*NonLine_Pk*M'+K*R0*eye(length(Z))*K'; % ����Э�������,Լɪ���ȶ���
           
           v_targetPostion=xk(1:3)';
           v_targetPostrue=xk(4:6)';
          %% ��¼���ε�������
          
          m_Alpha(k,:)=[v_targetPostion v_targetPostrue v_targetMagPolar];
        end
     end
    m_Alpha=m_Alpha(1:k,:);
    mean(m_Alpha)
    m_Alpha(end,:)
%% ����
    %% ��������
%      IterTimes=20;         % �������� 
%      Mu=4*10^13;                % ��������
%      AlphaTrue=[Target.Postion Target.Theta Target.Phi Target.Psi Target.MagPolar]; % ��ֵ
%      Alpha0=[0.3 0.2 -1 0 0 0 0.8 1 1];   % ������ֵ������Ϊx,y,z,�ȣ��գ��ף���x,��y,��z
%      m_Alpha=zeros(IterTimes+1,length(Alpha0));      % �洢��������
%      InvHz=zeros(size(Hz));                         % ���ݹ��������ɵ�ɨ�賡(��ǰֻ����z����)
%      HzGrad=zeros(IterTimes,length(Alpha0));        % ���ۺ������ݶ�
%      v_CostValue=zeros(IterTimes,1);                % ���ۺ����ڵ��������е�ֵ
%      m_Alpha(1,:)=Alpha0;
    

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
%      H=awgn(H,2000,'measured');
% 	fun = @(x)HzAlpha(x,Coils)-H';
%     options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt');
% [x,resnorm,residual,exitflag,output,lambda,jacobian]  = lsqnonlin(fun,Alpha0,[],[],options);
% resnorm
% output
% x
% AlphaTrue
% error=abs(AlphaTrue-x)
