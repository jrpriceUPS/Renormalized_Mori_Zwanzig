function [x,condition] = newton(f,grad_f,x0,tol)

err = tol + 1;
x = x0;
warning off MATLAB:nearlysingularMatrix
warning off MATLAB:singularMatrix
warning off MATLAB:illConditionedMatrix

while err > tol;
    y = f(x);
    grad_y = grad_f(x);

    change = grad_y\y;
    x = x - change;
    err = norm(change,2)/norm(x,2);
end

if sum(sum(isnan(grad_y)))+sum(sum(isinf(grad_y))) == 0
    condition = cond(grad_y);
else
    condition = NaN;
end