function [zfin, lam, delta] = potlanczos(x, coneidx, rho, g, f, ATA, AT, A, scale, xorig, vstart)
% Compute minimum eigen-value using Lanczos iteration for projected Hessian
% matrix
%
%       H = P * [- (g * g') / f + ATA) + diag(d * f / rho)] * P
%
% The matrix-vector multiplication is specially processed

if nargin < 11
    vstart = [];
end % End if

rng(24);
[~, n] = size(A);
maxiter = 150;

ncone = length(coneidx);
nmin = ceil(0.1 * ncone);
[~, xspmin] = mink(xorig(coneidx), nmin);
xspsmall = find(xorig(coneidx) < min(f, 1e-03));
xsp = union(xspmin, xspsmall);
xsp = xsp + n - length(coneidx);
% find(x(coneidx) < 1e-02) + n - length(coneidx);

if isempty(vstart)
if scale
%     v = x;
%     v(coneidx) = 1.0;
    v = g;
    v(coneidx) = v(coneidx) .* x(coneidx);
else
    v = g;
end % End if
else
    v = vstart;
end % End if

% v = randn(n, 1);
V = zeros(n, maxiter + 1);
H = zeros(maxiter + 1, maxiter);
v = v / norm(v);
V(:, 1) = v;
tol = 1e-06;

for k = 1:maxiter
    
    w = -MXv(v, coneidx, AT, A, ATA, x, f, g, rho, scale, xsp);
    
    wold = w;
    
    if (k > 1)
        w = w - H(k, k-1) * V(:, k-1);
    end % End if
    
    alp = w' * V(:, k);
    w = w - alp * V(:, k);
    H(k, k) = alp;
    
    if (norm(w) <= 0.99 * norm(wold) || 1)
        s = (w' * V(:,1:k))';
        w = w - V(:,1:k) * s;
        H(1:k, k) = H(1:k, k) + s;
    end % End if
    
    nrm = norm(w);
    v = w / nrm;
    V(:, k+1) = v;
    H(k+1,k) = nrm;  H(k,k+1) = nrm;
    
    if (rem(k, 5) == 0 || k == maxiter)
        
        Hk = H(1:k,1:k); Hk = 0.5 * (Hk + Hk');
        [Y, D] = eig(Hk);
        eigH  = real(diag(D));
        [~,idx] = sort(eigH);
        res_est = abs(H(k+1, k) * Y(k, idx(k)));
        
        if (res_est <= 0.01) || (k == maxiter)
            lam = eigH(idx(k));
            lam2 = eigH(idx(k-1));
            z = V(:, 1:k) * Y(:, idx(k));
            z2 = V(:, 1:k) * Y(:, idx(k - 1));
            zz = -MXv(z, coneidx, AT, A, ATA, x, f, g, rho, scale, xsp);
            
            res = norm(zz - lam * z);
            zfin = zz;
            zz = -MXv(z2, coneidx, AT, A, ATA, x, f, g, rho, scale, xsp);
            
            res2 = norm(zz - lam * z2);
            tmp = lam - lam2 - res2;
            
            if tmp > 0
                beta = tmp;
            else
                beta = eps;
            end % End if
            
            delta = min(res, res^2 / beta);
            
            if delta <= tol || (lam > 0 && lam - abs(delta) > 0)
                fprintf("Lanzos iteration %d \n", k);
                break;
            end % End if
            
        end % End if
    end % End if
end % End if


end % End function