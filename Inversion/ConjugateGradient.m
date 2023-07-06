function [i, xk1] = ConjugateGradient(f, grad, x0, iterations, tol)
%ConjugateGradient 非线性共轭梯度法
%   此处显示详细说明
% xlog = [];
% ylog = [];
% xdlog = [];
xk = x0;
c2 = 0.4;
flag = 0;
cs1=0;cs2=0;cs3=0;cs4=0;cs5=0;
fk = f(xk);
gk = grad(x0);
pk = -gk;

xk1 = xk;
gk1 = gk;
pk1 = pk;
%alpha1 = 2.0;
old_fval = f(x0);
old_old_fval = old_fval + norm(grad(x0)) / 2;
gnorm = max(abs(gk));
old_fval_back = old_fval;
old_old_fval_back = old_old_fval;

function [alpha,xkp1,pkp1,gfkp1,gnorm] = polak_ribiere_powell_step(alpha)
    xkp1 = xk + alpha * pk;
    gfkp1 = grad(xkp1);
    yk = gfkp1 - gk;
    beta_k = max(0, (yk*gfkp1')/(gk*gk'));
    pkp1 = -gfkp1 + beta_k * pk;
    gnorm = max(abs(gfkp1));
end
function n = descent_condition(alpha, xkp1, fp1, gfkp1)
    [cs1,cs2,cs3,cs4,cs5] = polak_ribiere_powell_step(alpha);
    alphak=cs1;xk1=cs2;pk1=cs3;gk1=cs4;gnorm=cs5;
    if gnorm <= 1e-5
        n = true;
    else
        n = (pk1*gk1') <= -0.01*(gk1*gk1');
    end
end
for i=1:iterations
%     xlog = [xlog xk'];
%     ylog = [ylog f(xk)];
%     fval = f(xk);
    %disp(i);
    
    old_fval_back = old_fval;
    old_old_fval_back = old_old_fval;
    [alpha, old_fval, old_old_fval] = StepLength3(f, grad, xk, pk, gk, old_fval, old_old_fval, 1e-4, 0.4, 1e100, 1e-100, 1e-14);
%     if isnan(alpha)
%         [alpha, old_old_fval] = StepLength2(f, grad, xk, pk, grad(xk), old_fval, old_old_fval, 1e-4, 0.4, 1e100, nan, nan, 20);
    %alpha = StepLength(f, grad, xk, pk, c2);
%     fprintf("iter=%d, x=%f,%f, f(x)=%f, gkpk=%f, alpha=%f\n",i, xk(1), xk(2), ylog(end), (gk*pk')/(gk*gk'), alpha);
    if ~isnan(alpha)
        if ~descent_condition(alpha,1,1,1)
            alpha = nan;
        end
    end
    if isnan(alpha)
        [alpha, old_fval, old_old_fval] = StepLength2(f, grad, xk, pk, grad(xk), old_fval_back, old_old_fval_back, 1e-4, 0.4, 1e100, @descent_condition, 1, 20);
    end
    if isnan(alpha)
%         ylog = [ylog nan];
%         fval = nan;
        break
    end
%     xk1 = xk + alpha*pk;
%     gk1 = grad(xk1);
    
    %下面是计算下一次共轭方向。《numerical optimization》P121
%     %beta_k1 = (gk1*gk1') / (gk*gk');
%     beta_k1 = (gk1*(gk1-gk)')/(gk*gk');
%     beta_k1 = max(0, beta_k1);
%     %disp(beta_k1);
%     pk1 = -gk1 +beta_k1 * pk;
%     fprintf("iter=%d, gk1=%f,%f, pk1=%f,%f, gk1pk1=%f, gk1pk=%f\n",i, gk1(1), gk1(2), pk1(1),pk1(2), (gk1*pk1')/(gk1*gk1'),(gk1*pk')/(gk*gk'));
    
    if xk1 ~= cs2
        [alpha, xk1, pk1,gk1,~] = polak_ribiere_powell_step(alpha);
    end
%     xdlog = [xdlog norm(xk1' - xlog(:,end))];
    curr_grad = abs(grad(xk1));
    if curr_grad < tol
        break;
    end
%     fprintf("iter=%d, y=%f, gk1=%f,%f, pk1=%f,%f, gkpk=%e, alpha=%f\n ",i, old_fval_back, gk(3), gk(4), pk(3), pk(4), gk*pk',alpha);
    %对于共轭梯度法，记录一下上次的步长alpha，用作下次寻找步长时的初值。《numerical
    %optimization》P59，式3.60上面的一个式子
    %alpha1 = (alpha*(gk*pk'))/(gk1*pk1');
    xk = xk1;
    gk = gk1;
    pk = pk1;
    
end
% fprintf('\niter=%d, x=%f,%f, f(x)=%f\n',i, xk(1), xk(2), ylog(end));
end

