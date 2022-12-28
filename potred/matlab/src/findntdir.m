function [m] = findntdir(A, b, c, x, y, s, tau)
% Compute inexact newton's direction

[m, n] = size(A);
mu = (x' * s) / n;
sigma = 1;

rp = A * x - b * tau;
rd = -A' * y - s + c * tau;
pobj = c' * x;
dobj = b' * y;

kappa = abs(pobj - dobj);

XSe = sqrt(x .* s); % n
D = sqrt(s) ./ sqrt(x); % n
Dinv = sqrt(x) ./ sqrt(s);

ADinv = A * sparse(1:n, 1:n, Dinv);
M = [speye(n), ADinv';
    ADinv, sparse(m, m)];

Dinvc = c ./ D;
DinvXinvrmu1 = XSe - (mu * sigma) ./ XSe;
rhs1 = [-Dinvc; b]; % m + n
rhs2 = [rd ./ D + DinvXinvrmu1; rp];
aux = [Dinvc; b]; % m + n

[L, DD, p] = ldl(M, 'vector');
p = p'; LT = L'; pinv = dsdpInvPerm(p)';

d1 = L \ rhs1(p);
d1 = DD \ d1;
d1 = LT \ d1;
d1 = d1(pinv);

d2 = L \ rhs2(p);
d2 = DD \ d2;
d2 = LT \ d2;
d2 = d2(pinv);

dtau = aux' * d2 + dobj - pobj - sigma * mu / tau;
dtau = - dtau / (kappa / tau - aux' * d1);

dxdy = d1 * dtau - d2;
dx = dxdy(1:n) ./ D;
dy = - dxdy(n+1:end);
dkappa = - kappa * dtau / tau - kappa + sigma * mu / tau;
ds = -D.* dxdy(1:n) - s + (mu * sigma) ./ x;

m = [dy; dx; ds; dtau];

end % End function