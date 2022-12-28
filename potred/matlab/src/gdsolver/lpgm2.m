function [lpsol] = lpgm2(A, b, ydim, maxiter)

[~, n] = size(A);
x = zeros(n, 1);
coneidx = ydim + 1:n;
x(coneidx) = 1.0;
aa = 1.0 / norz(A, 'fro')^2;

for i = 1:maxiter
    Axb = A * x - b;
    f = 0.5 * norm(Axb, 2)^2;
    g = A' * Axb;
    alpha = aa;
    x = x - alpha * g;
    x(coneidx) = max(x(coneidx), 0.0);
    if mod(i, 8) == 1
        aa = aa * 0.9;
    end % End if
    fprintf("%d  %5.3e \n", i, f);
end % End for

lpsol = x;

end % End function