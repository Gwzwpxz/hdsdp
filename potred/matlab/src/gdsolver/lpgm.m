function [lpsol, fvals] = lpgm(A, ydim, maxiter)

[~, n] = size(A);
x = zeros(n, 1);
coneidx = ydim + 1:n;
ncone = length(coneidx);
x(coneidx) = 1.0 / ncone;
lam = eigs(A' * A, 1, 'largestabs');
aa = 0.9 / lam;
fvals = zeros(maxiter, 1);
mu = 1.0;

for i = 1:maxiter
    
    Ax = A * x;
    f = 0.5 * norm(Ax, 2)^2;
    g = A' * Ax + mu * x;
    alpha = aa;
    x = x - alpha * g;
    x(coneidx) = projsplx(x(coneidx));
    fprintf("%d  %5.3e \n", i, f);
    fvals(i) = f;
    
    if mod(i, 32) == 1
        mu = mu * 0.1;
    end % End if
    
end % End for
lpsol = x;

end % End function