function [X,Xs,HzGrad,v_CostValue] = BFGS(Alpha,H,Coil,Plane,mu,precision,N)
%BFGS BFGS������С���˽�
%   ��������������ָ�����ȣ���BFGS�㷨�����С���˽⣬
%       �����
%           X��������ֵ
%           Xs,��������
%       ���룺    
%           Alpha��������ֵ
%           Coil,��Ȧ
%           H���۲�����(Ĭ��ȡHz��)
%           Plane���۲�ƽ��
%           mu,�������ӣ�������1��û�����Ż��������ӣ�
%           precision������
%           N,����������


%% ��ʼ������
X=Alpha(:);                           % x0��������ֵ
Xs=zeros(N,length(Alpha));            % N*m���󣬼�¼��������
D0=eye(length(Alpha));                % D0,��������hessian�������
v_CostValue=zeros(N,1);               % ���ۺ����ڵ��������е�ֵ
HzGrad=zeros(N,length(Alpha));        % ���ۺ������ݶ�
InvHz=zeros(size(H));                 % ���ݹ��������ɵ�ɨ�賡(��ǰֻ����z����)
%% ����gk
k=1;
        for i=1:length(Plane.v_y)          
            for j=1:length(Plane.v_x)
                 x=Plane.v_x(j);            % ��ǰɨ��ƽ���x����
                 y=Plane.v_y(i);            % ��ǰɨ��ƽ���y����
                 Coil.Postion(1:2)=[x y ];  % ������Ȧλ�ã�Ĭ�ϸ߶�Ϊ0���ֲ���  
                 [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,Alpha); % ��H��ǰλ�ó��ķֲ�
                 [temp1,temp2,Grad]=HGradAlpha(Coil,Alpha);   % ��ǰλ�õ��ݶ�
                 e=abs(InvHz(i,j))-abs(H(i,j));
                 HzGrad(k,:)=HzGrad(k,:)+2*e*Grad;
                 v_CostValue(k)=v_CostValue(k)+e^2;
            end
        end
        
gk=HzGrad(k,:)';
deltax=-mu*D0*gk;              % ţ�ٷ���
dis=sqrt(dot(deltax,deltax));
Xs(1,:)=X';
k=1;
while dis>precision &&k<N  
    k=k+1;
    X=X+deltax;
    Xs(k,:)=X';
        for i=1:length(Plane.v_y)          
            for j=1:length(Plane.v_x)
                 x=Plane.v_x(j);            % ��ǰɨ��ƽ���x����
                 y=Plane.v_y(i);            % ��ǰɨ��ƽ���y����
                 Coil.Postion(1:2)=[x y ];  % ������Ȧλ�ã�Ĭ�ϸ߶�Ϊ0���ֲ���  
                 [temp1,temp2,InvHz(i,j)]=HAlpha(Coil,X'); % ��H��ǰλ�ó��ķֲ�
                 [temp1,temp2,Grad]=HGradAlpha(Coil,X');   % ��ǰλ�õ��ݶ�
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

