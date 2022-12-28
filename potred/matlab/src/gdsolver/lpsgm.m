function [lpsol, fvals] = lpsgm(A, ydim, maxiter, p)

if nargin < 4
    p = 1;
    fprintf("Using 1 norm. \n");
end % End if

[~, n] = size(A);
x = zeros(n, 1);
coneidx = ydim + 1:n;
ncone = length(coneidx);
x(coneidx) = 1.0 / ncone;

nbuffer = 2;
Xbuffer = zeros(n, nbuffer);
econe = (1:nbuffer)' / ((nbuffer + 1) * nbuffer / 2);
rcount = 0;
aa = 1.0;

if p == 2
    aa = 1 / norm(A, 'fro');
end % End if

fvals = zeros(maxiter, 1);

xbest = x;
fbest = inf;

fthresh = norm(A * x, p);

for i = 1:maxiter
    
%     rcount = rcount + 1;
%     Xbuffer(:, rcount) = x;
    if rcount == nbuffer && 0
%         AX = A * Xbuffer;
%         cvx_begin quiet
%         cvx_solver gurobi
%         variable ads(nbuffer, 1)
%         minimize( norm(AX * ads, 1) )
%         subject to
%         Xbuffer(coneidx, :) * ads >= 0;
%         sum(ads) == 1.0;
%         cvx_end
%         ads = adsqp(AX, Xbuffer(coneidx, :), f * 0.1);    
        ads = econe;
        x = Xbuffer * ads;
        rcount = 0;
%         fprintf("%d  %5.3e \n", i, f);
%         aa = aa * 0.95;
    else
        Ax = A * x;
        
        if p == 1
            f = norm(Ax, 1);
            g = A' * sign(Ax);
            alpha = aa * f / norm(g)^2;
        else 
            f = norm(Ax, 1);
            g = A' * Ax / f;            
            alpha = aa;
        end % End if
        
        x = x - alpha * g;
        x(coneidx) = projsplx(x(coneidx));
%         fprintf("%d  %5.3e \n", i, f);
        
        if f < 0.95 * fthresh
            fthresh = f;
            aa = aa * 0.99;
        end % End if
        
%         if mod(i, 64) == 1
%             aa = aa * 0.5;
%         end % End if
        fvals(i) = f;
        
        if f < fbest
            fbest = f;
            xbest = x;
        end % End if
        
    end % End if
    
end % End for
lpsol = xbest;

end % End function