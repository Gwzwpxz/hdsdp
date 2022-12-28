function [x] = potProjObj(A, ydim, maxiter, mlp, nlp, Alp, blp, bscal, clp, cscal, E)
% Implement potential reduction for LPs
% The first ydim columns in A refer to y
% Projection is only done on partial variables

% warning off;
[~, n] = size(A); 
m = ydim;
ncone = n - m; 
yidx = 1:ydim;
coneidx = ydim + 1 : n;
projidx = coneidx;
rho = 1.1 * (ncone + sqrt(ncone));

rng(24);
x_prev = zeros(n, 1);
x_prev(coneidx) = 1.0;
x_prev(end) = 1;

if 0
    C = [x_prev'];
    CCTinv = full(inv(C * C'));
else
    C = [A(end, :)];
%     C = [x_prev';
%     A(end, :)];
%     A = A(1:end-1, :);
    CCTinv = full(inv(C * C'));
    cTx = -C(end, coneidx) * x_prev(coneidx);
    x_prev(yidx) = cTx * C(end, yidx)' / norm(C(end, yidx))^2;
    x_prev(end) = 1;
end % End if

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

betamax = 1.0;
beta = 1.0;

nrmb = norm(blp, 1);
nrmc = norm(clp, 1);

tic;
fprintf("%5s  %8s  %8s  %8s|  %8s  %8s  %10s  %10s %10s %10s\n",...
    "iter", "pres", "dres", "cpl",  "fval", "pot", "alphag", "alpham", "beta", "potred");

x_cum = x_pres;
logstar = "";
nbuffer = 8;
freq = 10;
Xbuff = zeros(n, nbuffer);
rcount = 0;

for i = 1:maxiter
    
    % Anderson acceleration
    if mod(i, freq) == 0
        rcount = rcount + 1;
        Xbuff(:, mod(rcount, nbuffer) + 1) = x_pres;
    end % End if
     
    if rcount >= nbuffer
        AX = A * Xbuff;
        andalp = adsqp(AX, Xbuff(coneidx, :), 0.0 * f);
        xand = Xbuff * andalp;
%         xand(projidx) = ncone * xand(projidx) / sum(xand(projidx));
        
        fand = 0.5 * norm(A * xand)^2;
        if fand < f
            
            potold = 1e+10;
            x_prev = x_pres;
            x_pres = xand;
            x_cum = x_cum .* x_pres;
            
%             rho = rho * 1.1;
            A = A * diag([ones(ydim, 1); x_pres(coneidx)]);
            ATA = A' * A;
            if size(C, 1) == 2
                C(1, coneidx) = C(1, coneidx) .* x_pres(coneidx)';
            end % End if
            C(end, coneidx) = C(end, coneidx) .* x_pres(coneidx)';
            CCTinv = full(inv(C * C'));
            x_prev(coneidx) = x_prev(coneidx) ./ x_pres(coneidx);
            x_pres(coneidx) = 1.0;
        end % End if
        
        rcount = 0;
    end % End if
    
    % Start potential reduction
    if recompute
        
        [f, g, ~] = fpot(A, ATA, x_pres);
        
        % Prepare momentum
        mk = x_pres - x_prev;
        
        % Gradient projection
        gk = g;
        gk(coneidx) = gk(coneidx) - (f ./ x_pres(coneidx)) / rho;
        % gk = (rho / f) * g;
        % gk(coneidx) = gk(coneidx) - (1 ./ x_pres(coneidx));
        
        gk = gk - C' * (CCTinv * (C * gk));
%         mk = zeros(n, 1);
        mk = mk - C' * (CCTinv * (C * mk));
        
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
    
%     if beta > 1
%         step = 1 / abs(min(d(coneidx) ./ x_pres(coneidx)));
%         step = min(0.9995 * step, 1.0);
%     else
%         step = 1.0;
%     end % End if
    
    xtmp = x_pres + d;
    
    if min(xtmp(coneidx)) < 0.0
        beta = beta * 0.25;
        recompute = false;
        continue;
    end % End if
    
    [ftmp, ~, ~] = fpot(A, ATA, xtmp);
    potnew = rho * log(ftmp) - sum(log(xtmp(coneidx)));
    potred = potnew - potold;    
    
%     aggstep = - 0.9995 / min(d(coneidx) ./ x_pres(coneidx));
%     decay = 0.8;
%     
%     while aggstep > 1.0
%         xtmp = x_pres + aggstep * d;
%         [linereduce, pottmp] = getpotreduce(rho, A, ATA, xtmp, potold, coneidx);
%         if linereduce < potred
%             potred = linereduce;
%             potnew = pottmp;
%             d = aggstep * d;
%             break;
%         end % End if
%         aggstep = aggstep * decay;
%     end % End while
    
    ratio = potred / mval;
    
    if isnan(ratio)
        ratio = 0;
    end % End if
    
    if f < 1e-20 || beta < 1e-05
        break;
    end % End if
    
    if (ratio < 0.0 || potred > 0)
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
    
    % Orthogonal projection to improve y
%     z1 = P * (-Q * x_pres(coneidx));
%     z2 = -C(end, coneidx) * x_pres(coneidx);
%     btmp = C(end, yidx)';
%     PPTinvb = PPT \ btmp;
%     PPTinvz1 = PPT \ z1;
%     nu = (btmp' * PPTinvz1 - z2) / (btmp' * PPTinvb);
%     x_pres(1:mlp) = PPTinvz1 - PPTinvb * nu;
%     x_pres(1:mlp) = PPT \ P * (-Q * x_pres(m+1:end));
    
    sol = x_pres .* E;
    sol(coneidx) = sol(coneidx) .* x_cum(coneidx);
    tau = sol(end);
    x = sol(mlp + 1: mlp + nlp) / (tau / bscal);
    y = sol(1:mlp) / (tau / cscal);
    s = sol(mlp + nlp + 1 : mlp + 2 * nlp) / (tau / cscal);
    
    pres = Alp * x - blp;
    dres = Alp' * y + s - clp;
    pobj = clp' * x;
    dobj = blp' * y;
    cpl = dobj - pobj;
    
    if mod(i, 500)
        fprintf("%5d | %8.2e %8.2e %8.2e | %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e  %+8.2e %1.1s %3.3f %+3.3e\n",...
            i, norm(pres) / (1 + nrmb), norm(dres) /  (1 + nrmc), abs(cpl) / (abs(pobj) + abs(dobj) + 1), f, potnew, alpha(1), alpha(2), beta, potred, logstar, toc, min(x_pres(coneidx)));
    end % End if
    
    logstar = "";
    
end % End for

x_pres(coneidx) = x_pres(coneidx) .* x_cum(coneidx);
x = x_pres;

end % End function