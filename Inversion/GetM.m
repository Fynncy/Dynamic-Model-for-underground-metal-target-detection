function [v_M] = GetM(Scene,r)
%�����������ݣ�����̽��Ŀ���ʼ״̬
%   ���: 
%       v_M���ż�������������Ԫ��������[M11,M22,M33,M12,M13,M23] ����ο�����
%   ����: 
%       Scene���������������в����Լ����շ�������
%       Scene.result.BxyzΪ���������ݣ�Scene.result.Bxyz0Ϊԭʼ���ݡ�
%       r����ǰλ��ʸ��
%       
%% ���ƴż���������M
Scene.model.metal.M;
% r=Scene.model.metal.postion(:); % ��ʼֵΪĿ����ʵλ��
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
end

