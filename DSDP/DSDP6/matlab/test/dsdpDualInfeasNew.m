function [S, y, kappa, tau, mu, pObj, reason, iter] = dsdpDualInfeasNew(A, b, C, y, S, Rd, mu, kappa, tau, dsdpParam)
% Implementation of dual scaling algorithm phase 1: dual infeasibility
% elimination

[n, ~] = size(S);
m = length(b);
maxiter      = dsdpParam{1};
tol          = dsdpParam{2};
nrmC         = norm(C, 'fro');
alphaphase1  = dsdpParam{14};
initstrategy = dsdpParam{7};
stepstrategy = dsdpParam{15};
ndash        = dsdpParam{20};
pweight      = 1.0; % dsdpParam{29};
prelax       = true;

pObj = inf;
step = 0;
delta = inf;
ub = 1e+07;
lb = -ub;
sl = y - lb * tau;
su = ub * tau - y;
dbigM = 1e+07;

nall = m * 2 + n;
binf = norm(b, 'inf');
bnorm1 = norm(b, 1);

mu = (ub * bnorm1 - dbigM * Rd(1, 1)) / (3 * nall);

rho = nall * 3;

for i = 1:maxiter
    
    nrmRd = sqrt(n) * abs(full(Rd(1, 1))) / (tau * (1 + nrmC));
    if nrmRd < tol / 10 && mu < 1e-05
        if ~ isinf(pObj)
            reason = "DSDP_PRIMAL_DUAL_FEASIBLE";
        else
            reason = "DSDP_PRIMAL_UNKNOWN_DUAL_FEASIBLE";
        end % End if
        break;
    end % End if
    
    dObj = b' * y;
    
    [M, u, asinv, ~, ~, ~, csinv, csinvcsinv, asinvrysinv, csinvrysinv, rysinv] = ...
        dsdpgetSchur(A, S, C, Rd, initstrategy);
    
    % Primal relaxation
    M = M + diag(sl.^-2 + su.^-2);
    % u = u + 1e+07 * (sl.^-2 + su.^-2);
    asinv = asinv - sl.^-1 + su.^-1;
    
    d1_2_3 = M \ [b, asinv, asinvrysinv];
    
    d1 = d1_2_3(:, 1);
    d2 = d1_2_3(:, 2);
    d3 = d1_2_3(:, 3);
    
    % Primal objective estimation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [mualpha, pObj, pinfeas] = dsdpgetmualpha(asinv, b, csinv, binf / min(abs(asinv)), ub, 1e-05);
    muestimate = (pObj - dObj - Rd(1, 1) * dbigM) / nall;
    if isinf(mu)
        mu = mualpha;
    else
        if mualpha > mu
            fprintf("Here\n");
        end % End
        mu = min([mu * 0.95, muestimate, mualpha]);
        % mu = max(min(mu * 0.85, muestimate), mu / 3);
    end % End if 
    
    dy  = 1 / mu * d1 - d2 + d3;
    dS  = Rd + C - dsdpgetATy(A, dy);
    % dkappa = - kappa + (mu / tau) - (kappa / tau) * dtau;
    % Compute stepsize
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    step = dsdpgetStepsize(S, dS, kappa, 1.0, tau, 0.0, stepstrategy, alphaphase1);
    
    if prelax
        dsu = -dy;
        dsl =  dy;
        tmp = - 1 / min(dsl ./ sl);
        if tmp > 0
            step = min(step, tmp);
        end % End if
        tmp = - 1 / min(dsu ./ su);
        if tmp > 0
            step = min(step, tmp);
        end % End if
    end % End if
    
    step = min(alphaphase1 * step, 1.0);
    
    if step < 1e-03
        reason = "DSDP_SMALL_STEP";
        break;
    end % End if
    
    % Potential reduction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delta = 0.0;
    % step = min(step, 0.9);
    [y, S, step] = dsdpInfeasPotReduction(A, b, C, y, dy, S, step, rho, pObj, Rd, ub, dbigM);
    Rd = Rd * (1 - step);
    
    % Take Newton step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    y = y + step * dy;
    % tau = tau + step * dtau;
    % kappa = mu / tau;
%     Rd = Rd * (1 - step);
%     S = - dsdpgetATy(A, y) + C * tau - Rd;
%     sl = y - lb * tau;
%     su = ub * tau - y;
    
    % Corrector
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [y, S] = dinfeaspotrdc(A, b * tau, C * tau, y, S, Rd, M, d2 * tau, mu, 4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Logging
    dObj = b' * y;
    fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n",...
        i, mualpha * csinv, dObj, nrmRd, kappa / tau, mu, step, pinfeas);
    
end % End for


% Corrector at the end
[y, S] = dinfeaspotrdc(A, b, C * tau, y, S, sparse(n, n), M, d2 * tau, mu, 12);

iter = i;

% Logging
fprintf("%3d  %10.2e  %10.2e  %8.2e  %8.2e  %8.2e  %8.2e  %8.2e \n",...
    i, pObj, b' * y / tau, nrmRd, kappa / tau, mu, step, delta);
showdash(ndash);
if reason == "DSDP_DUAL_INFEASIBLE"
    fprintf("Phase 1 certificates dual infeasibility \n");
elseif reason == "DSDP_SMALL_STEP"
    fprintf("Phase 1 ends due to small stepsize \n");
elseif reason == "DSDP_DUAL_FEASIBLE"
    fprintf("Phase 1 finds a dual feasible solution \n");
end % End if

end % End function