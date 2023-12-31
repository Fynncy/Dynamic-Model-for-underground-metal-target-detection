function [i, xk] = LevenbergMarquardt(f, grad, x0, iterations, tol)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明

xk = x0;
lamda = 0.01;
updataJ = 1;
tol2 = ones(1,length(xk))*tol;
for i = 1:iterations
    [g,J] = grad(xk);
    if updataJ == 1 
       H = J'*J;
       if i == 1
           e = f(xk);
       end
    end
    H_lm = H + (lamda * eye(length(xk),length(xk)));
    dp = -H_lm\g';  %警告: 矩阵为奇异值、接近奇异值或缩放错误。结果可能不准确。RCOND = NaN。
%     dp=-pinv(H_lm,1e-8)*g'; %报错 出错 pinv (line 18) [U,S,V] = svd(A,'econ'); 错误使用 svd  SVD 的输入不能包含 NaN 或 Inf。
%     dp=-lsqminnorm(H_lm,g'); %没有警告和报错，但是结果和警告出现奇异值的差不多
    xk1 = xk + dp';
    e_lm = f(xk1);
    
    if e_lm < e
       lamda = lamda/10;
       xk = xk1;
       e = e_lm;
%        disp(e);
       updataJ = 1;
    else
        updataJ = 0;
        lamda = lamda * 10;
    end
    if abs(g') < tol2 
        break;
    end
end
% fprintf("\niter=%d\n",i);
end
% inv和pinv
% 1.对于方阵A，如果为非奇异方阵，则存在逆矩阵inv(A)
% 2.对于奇异矩阵或者非方阵，并不存在逆矩阵，但可以使用pinv(A)求其伪逆
% 若A为非奇异矩阵，请不要使用pinv求逆，虽然计算结果相同，即
% inv( A ) = pinv( A )
% 但pinv的计算复杂度较高。
% 
% \和pinv （参考https://ww2.mathworks.cn/help/matlab/ref/pinv.html）
% A\B和pinv(A)*B 都是要解Ax=B
% 但是计算结果并不一样，\的更快更稳定，因为使用了高斯消去法
% 另外，官网提示 可以将应用于向量 b 的大部分 pinv（比如在 pinv(A)*b 中）替换为 lsqminnorm(A,b)，以计算线性方程组的最小范数最小二乘解。lsqminnorm 通常比 pinv 更有效，而且还支持稀疏矩阵。

%% position=[0 -1 -0.5];%目标位置
% \结果
% 【lm算法结果】
% 迭代次数：1000 运行时间：132.203
% 迭代初值：0  0 -5  0  0  0  0  0  0
% 迭代结果：-1.803001     -3.566669     -3.613827       428.923      72.39867      26.23171     -114.8834     -55.43245      33.30961
% 真实值：0+0i                    -1+0i                  -0.5+0i      0.016559+0.0003765i     0.019896+0.00036565i     0.034094+0.00031947i   -0.0051674+1.6808e-05i   -0.0089501+2.9112e-05i     0.012295-3.9992e-05i
% 误差：1.803001      2.566669      3.113827      428.9064      72.37877      26.19762      114.8782       55.4235      33.29731
% 【目标函数值】
% 迭代初值处：87.022
% 收敛值处：48.5105
% 真实值处：0.51161
% 【定位误差/m】
% x,y,z=1.803      2.5667      3.1138
% 【主轴极化率误差/百分比】
% x,y,z=8.6964 47.289 471.4974
% 【姿态误差/角度】
% theta,phi,psi=51.5243 33.9129 0
% 
% lsqminnorm结果
% test_duodian
% 【lm算法结果】
% 迭代次数：1000 运行时间：140.405
% 迭代初值：0  0 -5  0  0  0  0  0  0
% 迭代结果：-1.805936     -3.574142     -3.613576      431.4847      72.78505      26.12863     -115.3621     -55.99763      33.75068
% 真实值：0+0i                    -1+0i                  -0.5+0i      0.016559+0.0003765i     0.019896+0.00036565i     0.034094+0.00031947i   -0.0051674+1.6808e-05i   -0.0089501+2.9112e-05i     0.012295-3.9992e-05i
% 误差：1.805936      2.574142      3.113576      431.4682      72.76515      26.09453      115.3569      55.98868      33.73839
% 【目标函数值】
% 迭代初值处：87.0221
% 收敛值处：48.5097
% 真实值处：0.49914
% 【定位误差/m】
% x,y,z=1.8059      2.5741      3.1136
% 【主轴极化率误差/百分比】
% x,y,z=8.3446 47.709 474.2743
% 【姿态误差/角度】
% theta,phi,psi=51.5288 33.794 0
% 
