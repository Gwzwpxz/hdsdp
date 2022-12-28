function [data] = presgdtrans(fpath)
% Transform LP into SGD standard form

prob = gurobi_read(fpath);
% prob = gurobi_presolve(prob);
eqidx = find(prob.sense == '=');
geqidx = find(prob.sense == '>');
leqidx = find(prob.sense == '<');

c = prob.obj;
lb = prob.lb;
ub = prob.ub;
lb(lb == -inf) = -1e+06;
ub = ub - lb;

A = prob.A;
[~, n] = size(A);
b = prob.rhs - A * lb;

boundid = find(ub ~= +inf);
nbound = length(boundid);
G = sparse(1:nbound, boundid, ones(nbound, 1), nbound, n);
data.A = A(eqidx, :);
data.b = b(eqidx);
data.G = [ A(geqidx, :);
          -A(leqidx, :);
          -G];
data.g = [b(geqidx); -b(leqidx); -ub(boundid)];
data.c = c;

end % End function
