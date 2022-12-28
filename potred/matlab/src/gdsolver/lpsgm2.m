function [lpsol] = lpsgm2(A, b, ydim, x0, maxiter, p)

[~, n] = size(A);
coneidx = ydim + 1:n;
if isempty(x0)
    x = zeros(n, 1);
    x(coneidx) = 1.0;
else
    x = x0;
end % End if

aa = 1.0;
fbest = inf;
xbest = x;

for i = 1:maxiter
    
    Axb = A * x - b;
    Axb(end) = Axb(end) * 1000;
    
    if p == 1
        f = norm(Axb, 1);
        g = A' * sign(Axb);
    else
        f = norm(Axb, 2);
        g = A' * Axb / f;
    end % End if
    
    alpha = aa * f / norm(g)^2;
    x = x - alpha * g;
    x(coneidx) = max(x(coneidx), 0.0);
    %     if mod(i, 8) == 1
    %         aa = aa * 0.9;
    %     end % End if
%     fprintf("%d  %5.3e \n", i, f);
    
    if f < fbest
        fbest = f;
        xbest = x;
    end % End if
end % End for

lpsol = xbest;

end % End function