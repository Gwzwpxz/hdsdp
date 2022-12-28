function [data] = presgdtrans2(fpath)
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
b = prob.rhs - A * lb;

data.objcon = c' * lb;
data.A = A(eqidx, :);
data.b = b(eqidx);
data.G = [ A(geqidx, :);
          -A(leqidx, :)];
data.g = [b(geqidx); -b(leqidx)];
data.c = c;
data.u = ub;

end % End function
