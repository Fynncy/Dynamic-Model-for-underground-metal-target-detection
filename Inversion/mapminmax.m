function [output] = mapminmax(v_data,l,r)
%MAPMINMAX 此处显示有关此函数的摘要
%   此处显示详细说明
mi=min(v_data(:));
mx=max(v_data(:));
output=(v_data-mi)/(mx-mi);
output=(output)*(r-l)+l;
end

