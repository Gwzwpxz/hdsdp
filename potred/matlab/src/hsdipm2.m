function [x, y, s, kappa, tau] = hsdipm2(A, b, c)
% Implement a primal-dual interior point method using HSD embedding

warning off;
[m, n] = size(A);
x = ones(n, 1);
y = zeros(m, 1);
s = ones(n, 1);
tau = 1;
kappa = 1;

sigma = 0.7;

cnrm = norm(c, 1) + 1;
bnrm = norm(b, 1) + 1;
fprintf("%4s %8s %8s %8s %8s %8s \n", "iter", "pobj", "dobj", "pinf", "dinf", "mu");

for i = 1:100
    
    mu = (x' * s + kappa * tau) / (n + 1);
    
    % Set up residuals
    rp = A * x - b * tau;
    rd = -A' * y - s + c * tau;
    pobj = c' * x;
    dobj = b' * y;
    
    pinf = norm(rp) / (tau * bnrm);
    dinf = norm(rd) / (tau * cnrm);
    
    if max([pinf, dinf, mu]) < 1e-08
        break;
    end % End if
    
    fprintf("%4d %8.1e %8.1e %8.1e %8.1e %8.1e \n", i, pobj / tau, dobj / tau, pinf, dinf, mu);
    
    xinvs = s ./ x;
    
    M = [sparse(1:n, 1:n, xinvs), A';
         A,                  sparse(m, m)];
    
    Xinvrmu1 = s - mu * sigma * x.^-1; 
    rhs1 = [-c; b]; % m + n
    rhs2 = [rd + Xinvrmu1; rp];
    aux = [c; b]; % m + n
    
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
    dx = dxdy(1:n);
    dy = - dxdy(n+1:end);
    dkappa = - kappa * dtau / tau - kappa + sigma * mu / tau;
    ds = -dxdy(1:n) .* xinvs - Xinvrmu1;
    
    alpha = 0.995 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau]));
    
    x = x + alpha * dx;
    y = y + alpha * dy;
    s = s + alpha * ds;
    kappa = kappa + alpha * dkappa;
    tau = tau + alpha * dtau;
    
    simp = sum(x) + sum(s) + kappa + tau;
    simp = simp / (2 * n + 2);
    x = x / simp;
    y = y / simp;
    s = s / simp;
    kappa = kappa / simp;
    tau = tau / simp;
    
end % End for

end % End function