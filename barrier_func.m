function [B,g,h] = barrier_func(funs, x, t)
[f, cons, con_grad, con_hess] = funs(x);
[~,k] = size(con);
[n,~] = size(con_grad{1});
g = zeros(n,1);
h = zeros(n,n);
if all(cons <= 0)
    %x is feasible
    B = t*f - sum( log(-cons));
    for i = 1:k
       g = g + (-1 / con(i))*con_grad{i};
       h = h + (1 / con(i)^2)*con_grad{i}*con_grad{i}' + (-1 / con(i))*con_hess{i};
    end
else
    B = inf;
    g = inf*ones(n,1);
    h = inf*ones(n,n);
end