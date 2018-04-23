function [x,f,iters,grad_norms] = newtMethBarrier(funs, t, x0, tol, maxit, alpha, beta)

[f0, g, H] = barrier_func(funs, x0, t);

iters = 1;
x = x0;
f = f0;
grad_norms = norm(g);
while((iters <= maxit) && norm(g) > tol)
    delta_x = -H \ g;        % a gradient step
    step_size = BTLS(funs, delta_x, x, alpha, beta);
    x = x + step_size*delta_x;
    [f, g, H] = fun(x);
    iters = iters + 1;
    grad_norms = [grad_norms, norm(g)];
end

