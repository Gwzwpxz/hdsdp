clear;
% clc;
% close all;

fname = fullfile('./data', 'p_afiro.mps');
% fname = '/Users/gaowenzhi/Desktop/mipheurs/l2cta3d/shift.mps.gz';

data = preprocess(fname);
rng(24);

A = data.A;
b = data.b;
c = data.c;

% savepotdata(data);

linesearch = false;
neweigs = false;

[m, n] = size(A);

p = 1; 
nrmb = norm(b, p); 
nrmc = norm(c, p);

bscal = 1; nrmb; 
cscal = 1; nrmc;
b = b / bscal; 
c = c / cscal;
% model.rhs = b;
% model.obj = c;
% grbsol = gurobi(model);

% HSDAA = [sparse(m, m), A, sparse(m, n), sparse(m, 1), -b;
%          -A',         sparse(n, n),  -speye(n), sparse(n, 1), c;
%          b',          -c',  sparse(1, n), -1, 0];
     
HSDAA = [sparse(m, m), A, sparse(m, n), -b;
         -A',         sparse(n, n),  -speye(n), c;
         b',          -c',  sparse(1, n), 0];
     
% HSDAA = HSDAA / 100;
[D, E, HSDA] = ruizscale(HSDAA, 0);
E = 1./ E;
D = 1./ D;
[D2, E2, HSDA] = pcscale(HSDA, 0);
D = D .* D2;
E = E .* E2;

pdratio = norm(c) / norm(b);
HSDA(1:m, :) = HSDA(1:m, :); 
HSDA(m+1:end, :) = HSDA(m+1:end, :); 
% HSDA(end, :) = HSDA(end, :) * 10; 
% lpsol = potreduceLp(HSDA, m, 5000, false, linesearch, neweigs, 1);
% lpsol = potProjObj(HSDA, m, 80000, m, n, data.A, data.b, bscal, data.c, cscal, E);
% [lpsol, fvals] = lpsgm(HSDA, m, 1000, 1);
lpsol = potRecur(HSDA, m, 10000, m, n, A, b, c, E);
% HSDA(end, :) = HSDA(end, :);
% lpsol = lpsgm2(HSDA(:, 1:end-1), -HSDA(:, end), m, x0, 1000, 2);
sol = lpsol .* E; 
tau = lpsol(end);

% tau = sol(end);
y = sol(1:m);
y = y / (tau / cscal); 
s = sol(m + n + 1 : m + 2 * n) / (tau / cscal);
x = sol(m + 1: m + n) / (tau / bscal);

pobj = data.c' * x;
dobj = data.b' * y;
pres = norm(data.A * x - data.b);
dres = norm(data.A' * y + s - data.c);

fprintf("%20s %10.3e, %10.3e  %10.3e  %10.3e \n", '', pobj, ...
    pres / (1 + norm(data.b, 1)), dres / (norm(data.c, 1) + 1), ...
    abs(pobj - dobj) / (abs(pobj) + abs(dobj) + 1));
