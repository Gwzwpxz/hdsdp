function [x, y, s, w] = sglpsolve2(c, A, b, G, g, u, x0, y0, maxit, quiet)
% Solve LP using penalty based method using a sharp objective
% min  c' * x
% s.t.  A * x  = b
%       G * x >= g
%      u >= x >= 0

[meq, n] = size(A);
[mineq, n2] = size(G);

if meq == 0 || mineq == 0
    if meq > 0
        G = sparse(0, n);
        n2 = n;
    else
        A = sparse(0, n2);
        n = n2;
    end % End if
end % End if
      
assert( n == n2 );
ubid = find(u ~= inf);
nub = length(ubid);
u = u(ubid);

ineqidx = meq+1:meq+mineq;
K = [A; G];
q = [b; g];

% [D, E, K] = pcscale(K, 1); 

m = meq + mineq;

if isempty(x0)
    x = ones(n, 1);
    y = zeros(m, 1);
    w = zeros(nub, 1);
else
    x = x0;
    y = y0;
end % End if

pr = 1;
dr = 1;
pdr = 1;
aa = 1;

fbest = inf;
xbest = x0;
ybest = y0;
freq = 10;

for i = 1:maxit
    
    % Subgradient
    % Primal x
    Kx = K * x - q;
    Kx(ineqidx) = min(Kx(ineqidx), 0);
    dx = K' * sign(Kx);
    
    % Dual y and w
    KTy = K' * y - c;
    KTy(ubid) = KTy(ubid) - w;
    KTy = max(KTy, 0);
    sKTy = sign(KTy);
    dy = K * sKTy;
    dw = -sKTy(ubid);
    
    % Compl.
    cpl = c' * x - q' * y + u' * w;
    dx = pr * dx + pdr * c * sign(cpl);
    dy = dr * dy - pdr * q * sign(cpl);
    dw = dr * dw + pdr * u * sign(cpl);
    
    % Residual
    pres = sum(abs(Kx));
    dres = sum(abs(KTy));
    cres = abs(cpl);
    fval = pr * pres + dr * dres + pdr * cres;
    
    if fval < 1e-08
        break;
    end % End if
    
    if fval < fbest
        fbest = fval;
        xbest = x;
        ybest = y;
        wbest = w;
    end % End if
    
    % Polyak
    nrmdx = norm(dx);
    nrmdy = norm(dy);
    nrmdw = norm(dw);
%     alpha = max(aa * fval, 1.5 * fbest) / (nrmdx^2 + nrmdy^2 + nrmdw^2);
    alpha = aa * fbest / (nrmdx^2 + nrmdy^2 + nrmdw^2);
%     alpha = aa * fval / (nrmdx^2 + nrmdy^2 + nrmdw^2);
    
    % Update
    x = x - alpha * dx;
    y = y - alpha * dy;
    w = w - alpha * dw;
    
    % Projection
    x = max(x, 0);
    x(ubid) = min(x(ubid), u);
    y(ineqidx) = max(y(ineqidx), 0);
    w = max(w, 0);
    
    if 1 % mod(i, 100) == 1 && ~quiet
        fprintf("%4d %5.2e %5.2e %5.2e | %5.2e \n", i, pres, dres, cres, fval);
    end % End if

    if mod(i, freq) == 1
        aa = aa * 0.95;
    end % End if

end % End for

x = xbest;
y = ybest;
ww = zeros(n, 1);
ww(ubid) = wbest;
w = ww;
s = max(c - K' * y + w, 0);

end % End function