clear; clc; close all;
m = 1000;
n = 10000;

A = rand(m, n);
k = 64;
nc = 90;

X = randn(n, k);
AX = A * X;
alpha = adsqp(AX, X(n - nc + 1:end, :));

model.Q = 0.5 * sparse(AX' * AX);
model.A = [X(n - nc + 1:end, :);
           ones(1, k)];
model.A = sparse(model.A);
model.rhs = zeros(nc + 1, 1);
model.sense = '>';
model.lb = -inf(k, 1);
model.rhs(end) = 1.0;

sol = gurobi(model);
