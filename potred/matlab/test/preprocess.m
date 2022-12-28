function [data] = preprocess(fpath)
% generate params
% prob = mpsread(fpath);

% wrap up to the format for abip
% Aeq = prob.Aeq;
% Aineq = prob.Aineq;
% beq = prob.beq;
% bineq = prob.bineq;

prob = gurobi_read(fpath);
% prob = gurobi_presolve(prob);
eqidx = find(prob.sense == '=');
geqidx = find(prob.sense == '>');
leqidx = find(prob.sense == '<');

Aeq = prob.A(eqidx, :);
beq = prob.rhs(eqidx, :);
Aineq = [prob.A(leqidx, :);
         -prob.A(geqidx, :)];
bineq = [prob.rhs(leqidx); -prob.rhs(geqidx)];

% translate to equality constraints
% length of constraints
[m2, ~] = size(Aineq);
[m1, n] = size(Aeq);

% inspect bounds
% if unbounded below set -1e6
lb = (prob.lb > - inf) .* prob.lb;
lb(isnan(lb)) = -1e+06;
lb = lb + (prob.lb == - inf) .* (-1e+06);

% consider x <= ub
ub = prob.ub;
% record non-inf ub
idxub = ub < inf;
m3 = sum(idxub);
D = spdiags(idxub(:), 0, n, n);
D = D(idxub, :);
% produce a bound induced rhs
brhs = prob.ub(idxub) - lb(idxub); 

% LHS matrix
A = [Aeq, sparse(m1, m2 + m3);
    Aineq, speye(m2), sparse(m2, m3);
    D, sparse(m3, m2), speye(m3)];
% RHS
b = [beq - Aeq * lb;
    bineq - Aineq * lb;
    brhs];
c = [prob.obj; sparse(m2 + m3, 1)];
% ub is not needed now.
% set 0s as lower bound lb
lb = full(sparse(n + m2 + m3, 1));
% some statistic of A

data.A = sparse(A);
data.b = full(b);
data.c = full(c);
data.lb = lb;
data.objcon = prob.obj' * prob.lb;

