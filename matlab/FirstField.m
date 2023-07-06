function [Hxyz] = FirstField(I,R,frequency,postion)
%������Ȧ����������ռ���ĳ��Ĵų�ʸ��
%   ���: 
%       Hxyz,3x1��������ΪHx,Hy,Hz
%   ����: 
%       I������
%       R����Ȧ�뾶
%       frequency��Ƶ��
%       postion���۲��λ��,3x1����������Ϊx,y,z
A=R^2*I;
k=2*pi*frequency/(3*10^8);
x=postion(1);
y=postion(2);
z=postion(3);
r=sqrt(sum(postion.^2));
r_xy=sqrt(x^2+y^2); %wz19.12.19

 sint=r_xy/r;
 cost=z/r;
 sinf=y/r_xy;
 cosf=x/r_xy;
 r3=r^3;
 r2=r^2; 
 He=exp(-1j*k*r);
        
 Hr=A*cost*0.5*(1/r3+1j*k/r2);
 Hr=Hr*He;
 Ht=A*sint*0.25*((1/r3-(k^2)/r)+1j*k/r2);
 Ht=Ht*He;
 if r_xy<realmin  % Ŀ�����z����
     Hxyz(1)=0;
     Hxyz(2)=0;
 else 
      Hxyz(1)=Hr*sint*cosf+Ht*cost*cosf;
      Hxyz(2)=Hr*sint*sinf+Ht*cost*sinf;
     
 end
 %Hxyz(3)=Hr*cost+Ht*sint;
 Hxyz(3)=Hr*cost-Ht*sint; %wz19.12.19
% Hxyz=abs(Hxyz);
end

