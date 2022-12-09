clear; clc; close all;

load('ads.mat');
AX = data.AX;
X = data.X;

k = 16;
nc = 58;

alpha = adsqp(AX, X(19:end, :));

model.Q = 0.5 * sparse(AX' * AX);
model.A = [X(19:end, :);
           ones(1, k)];
model.A = sparse(model.A);
model.rhs = zeros(nc + 1, 1);
model.sense = '>';
model.lb = -inf(k, 1);
model.rhs(end) = 1.0;

sol = gurobi(model);
