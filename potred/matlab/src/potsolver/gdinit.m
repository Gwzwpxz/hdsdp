function [u] = gdinit(M, v, m, mu, maxiter)
% Initialize the potential solver
[~, n] = size(M);
u = zeros(n, 1);
coneidx = m + 1:n;
u(coneidx) = 1.0;
aa = 1.0 / svds(M, 1, 'largest')^2;

for i = 1:maxiter
    Mxv = M * u - v;
    f = 0.5 * norm(Mxv, 2)^2 + 0.5 * mu * norm(u)^2;
    g = M' * Mxv + mu * u;
    u = u - aa * g;
    u(coneidx) = max(u(coneidx), 0.0);
    
    nrmg = norm(g);
%     fprintf("%d  %5.3e %5.3e\n", i, f, nrmg);
    
    if nrmg < 1e-05
        break;
    end % End if
    
end % End for

end % End function