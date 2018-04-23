function [f,con, con_grad, con_hess] = test_problem(x)

load('quadratics')
[~, k] = size(A);
 
con = zeros(1,k);
con_grad = cell(1,k);
con_hess = cell(1,k);

f = sum(x);

for i = 1:k
   con(i) = x'*A{k}*x;
   con_grad{i} = A{k}*x + A{k}'*x;
   con_hess{i} = A{k} + A{k}';
end

end