function [x] = potRecur(A, ydim, maxiter, mlp, nlp, Alp, blp, clp, E)
% Implement potential reduction for LPs
% The first ydim columns in A refer to y

% warning off;
[~, n] = size(A); m = ydim;
ncone = n - m; coneidx = ydim + 1 : n;
projidx = coneidx;
% projidx = n - 1 : n;
nproj = length(projidx);
rho = 1.1 * (ncone + sqrt(ncone));
e = ones(nproj, 1);

rng(24);
x_prev = zeros(n, 1);
x_prev(coneidx) = 1.0;

AT = A';
ATA = A' * A;
% One step potential reduction
[f, ~] = fpot(A, ATA, x_prev);
potold = rho * log(f) - sum(log(x_prev(coneidx)));
x_pres = x_prev;
recompute = true;

PT = A(mlp + 1:end, 1:mlp);
PPT = PT' * PT;
P = PT';
Q = A(mlp + 1:end, mlp+1:end);

fvals = zeros(maxiter, 1);
potrds = zeros(maxiter, 1);

betamax = 1.0;
beta = betamax;

tic;
fprintf("%5s  %8s  %8s  %8s|  %8s  %8s  %10s  %10s %10s %10s\n",...
    "iter", "pres", "dres", "cpl",  "fval", "pot", "alphag", "alpham", "beta", "potred");

x_cum = x_pres;

logstar = "";
allowcurv = 0;
usecurvature = 1;
curvinterval = 0;
ncurvs = 0;

nbuffer = 32;
freq = 2000000;
Xbuff = zeros(n, nbuffer);
rcount = 0;

for i = 1:maxiter
    
    % Anderson acceleration
    if mod(i, freq) == 0
        rcount = rcount + 1;
        Xbuff(:, mod(rcount, nbuffer) + 1) = x_pres;
    end % End if
     
    if rcount >= nbuffer
%         Xbuff(:, 1) = randn(n, 1);
        AX = A * Xbuff;
        andalp = adsqp(AX, Xbuff(coneidx, :), 0.0 * f);
        xand = Xbuff * andalp;
        xand(projidx) = ncone * xand(projidx) / sum(xand(projidx));
        
        fand = 0.5 * norm(A * xand)^2;
        if fand < f
            
            potold = 1e+10;
            x_prev = x_pres;
            x_pres = xand;
            x_cum = x_cum .* x_pres;
            
            A = A * diag([ones(ydim, 1); x_pres(coneidx)]);
            AT = A';
            ATA = A' * A;
            
%             PT = A(mlp + 1:end, 1:mlp);
%             PPT = PT' * PT;
%             P = PT';
%             Q = A(mlp + 1:end, mlp+1:end);
            
            x_prev(coneidx) = x_prev(coneidx) ./ x_pres(coneidx);
            x_pres(coneidx) = 1.0;
        end % End if
        
        rcount = 0;
    end % End if
    
    % Start potential reduction
    if recompute
        
        [f, g, ~] = fpot(A, ATA, x_pres);
        fvals(i) = f;
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        if i >= 100
            mk = mnt;
        end % End if
        
        % Gradient projection
        gk = g;
        gk(coneidx) = gk(coneidx) - (f ./ x_pres(coneidx)) / rho;
        
        if usecurvature && allowcurv
            logstar = "*";
            method = "direct";
            xorig = x_pres;
            xorig(coneidx) = xorig(coneidx) .* x_cum(coneidx);
            [mk, ~] = findnegacurv(x_pres, m, coneidx, projidx, rho, g, f, ATA, AT, A, [], method, xorig);
            usecurvature = false;
            ncurvs = ncurvs + 1;
        end % End if
        
        gk(projidx) = gk(projidx) - e * sum(gk(projidx)) / nproj;
        mk(projidx) = mk(projidx) - e * sum(mk(projidx)) / nproj;
        
        % Prepare hessian
        nrmgk = norm(gk);
        nrmmk = norm(mk);
        
        if nrmmk > 0.0
            mk = mk / nrmmk;
        end % End if
        
        gk = gk / nrmgk;
        nrmgk = nrmgk * rho / f;
        Agk = A * gk; Amk = A * mk;
        
        gTgk = g' * gk; gTmk = g' * mk;
        xinvgk = gk ./ x_pres; xinvmk = mk ./ x_pres;
        xinvgk(1:m) = 0; xinvmk(1:m) = 0;
        
        gkXXgk = norm(xinvgk)^2;
        mkXXmk = norm(xinvmk)^2;
        gkXXmk = xinvgk' * xinvmk;
        
        gkHgk = -rho * (gTgk / f)^2 + gkXXgk + rho * norm(Agk)^2 / f;
        mkHmk = -rho * (gTmk / f)^2 + mkXXmk + rho * norm(Amk)^2 / f;
        mkHgk = -rho * (gTmk * gTgk) / f^2 + gkXXmk + rho * Amk' * Agk / f;
        
        gkTgk = gk' * gk * nrmgk; gkTmk = gk' * mk * nrmgk;
        
        H = [gkHgk, mkHgk;
             mkHgk, mkHmk];
        h = [gkTgk; gkTmk];
        M = [gkXXgk, gkXXmk;
             gkXXmk, mkXXmk];
        
    end % End if
    
    [alpha, mval] = subtrust(H, h, M, beta^2 / 2.5, 1e-10);
    d = alpha(1) * gk + alpha(2) * mk;
    
    xtmp = x_pres + d;
    [ftmp, ~, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;
    
    potrds(i) = potred;
    
    ratio = potred / mval;
    
    if isnan(ratio)
        ratio = 0;
    end % End if
    
    if f < 1e-20 || beta < 1e-05
        break;
    end % End if
    
    if (ratio < 0.2 || potred > 0)
        beta = beta * 0.25;
        recompute = false;
        continue;
    elseif ratio > 0.75
        beta = min(beta * 2, betamax);
    else
        
    end % End if
    
    recompute = true;
    x_prev = x_pres;
    x_pres = x_pres + d;
    potold = potnew;
    
%     x_pres(1:mlp) = PPT \ P * (-Q * x_pres(m+1:end));
    
    curvinterval = curvinterval + 1;
    if potred > -50 && curvinterval > 100
        usecurvature = true;
        curvinterval = 0;
    end % End if
    
    sol = x_pres .* E;
    sol(coneidx) = sol(coneidx) .* x_cum(coneidx);
    tau = sol(end);
    x = sol(mlp + 1: mlp + nlp) / tau;
    y = sol(1:mlp) / tau;
    s = sol(mlp + nlp + 1 : mlp + 2 * nlp) / tau;
    
    pres = Alp * x - blp;
    dres = Alp' * y + s - clp;
    cpl = blp' * y - clp' * x;
    
    mnt = findntdir(Alp, blp, clp, x, y, s, 1);
    mnt = mnt ./ E;
    
    if mod(i, 500)
        fprintf("%5d | %8.2e %8.2e %8.2e | %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f %+3.3e\n",...
            i, norm(pres), norm(dres), norm(cpl), f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc, min(x_pres(coneidx)));
    end % End if
    
    logstar = "";
    
end % End for

fprintf("Curvature is computed %d times. \n", ncurvs);
x_pres(coneidx) = x_pres(coneidx) .* x_cum(coneidx);
x = x_pres;
% semilogy(fvals, 'LineWidth', 3);
% hold on;

end % End function