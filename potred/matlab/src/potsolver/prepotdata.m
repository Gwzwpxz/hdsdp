function [data] = prepotdata(fpath)
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
G = [-A(geqidx, :);
     A(leqidx, :)];
 
data.objcon = c' * lb;
data.A = [A(eqidx, :), sparse(length(eqidx), size(G, 1));
          G, speye(size(G, 1))];
data.b = [b(eqidx); -b(geqidx); b(leqidx)];
data.c = [c; zeros(size(G, 1), 1)];
data.u = [ub; +inf(size(G, 1), 1)];

end % End function
