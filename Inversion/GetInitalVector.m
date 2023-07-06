function [v_inital] = GetInitalVector(Scene)
%�����������ݣ�����̽��Ŀ���ʼ״̬
%   ���: 
%       v_inital����ʼ״̬��[x,y,z,M11,M22,M33,M12,M13,M23] ����ο�����
%   ����: 
%       Scene���������������в����Լ����շ�������
%       Scene.result.BxyzΪ���������ݣ�Scene.result.Bxyz0Ϊԭʼ���ݡ�
%       
v_inital=zeros(9,1);
%% ����ˮƽλ��
x=sum(sum(Scene.result.Bxyz.*Scene.dataconf.m_x))/sum(sum(Scene.result.Bxyz));
y=sum(sum(Scene.result.Bxyz.*Scene.dataconf.m_y))/sum(sum(Scene.result.Bxyz));

v_inital(1:2)=[x y]';
%% �������
area=Scene.dataconf.linespace*Scene.dataconf.interval;
z=-2*sqrt(FWHM(Scene.result.Bxyz,area)/pi);
v_inital(3)=z;
%% ���ƴż���������M
Scene.model.metal.M;
r=[x y z]'; %%��ʼλ��ʸ���ɿռ���Ӧһ�׾غͰ�岨��ȷ��
% r=Scene.model.metal.postion(:); % ��ʼֵΪĿ����ʵλ��
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
              Pk=[Scene.dataconf.v_x(j) Scene.dataconf.v_y(i) Scene.dataconf.height]';  %��ǰ��Ȧλ��
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

