function [output] = mapminmax(v_data,l,r)
%MAPMINMAX �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
mi=min(v_data(:));
mx=max(v_data(:));
output=(v_data-mi)/(mx-mi);
output=(output)*(r-l)+l;
end

