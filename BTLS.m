function step_size = BTLS(fun, delta_x, x, alpha, beta, t)
%code for backtracking line search.  The step size, t, is initialized at
%one.  The search direction, delta_x, is passed in as an argument.
step_size = 1;
[f0, g0, ~] = barrier_func(fun, x, t);
[f1, ~, ~] = barrier_func(fun, x + step_size*delta_x, t);
f_star = f0 + alpha*step_size*g0'*delta_x;
while(f1 > f_star)
    step_size = beta*step_size;
    [f1, ~, ~] = barrier_func(funs, x + step_size*delta_x, t);
    f_star = f0 + alpha*step_size*g0'*delta_x;
end