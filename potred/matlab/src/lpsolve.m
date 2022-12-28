function [x, y, s, w] = lpsolve(A, b, c, ub, objcon, params)
% Solve an LP using first order potential reduction

% Get problem statistics
[m, n] = size(A);


% Prepare scaling and cone location
[scal, M, v, ubidinall, ubsub, coneid] = ...
    lpscale(A, b, c, ub, params.ruiz, params.l2scale, params.method);
scal.objcon = objcon;

% Initialize iterations and balance residuals
[u, z] = initpotlp(m, n, M, v, ubidinall, ubsub, coneid, scal);
[~, rp, rd, ~] = potresi(m, M, u, v, 1.0);
[~, ~, presrel, dresrel, ~, gaprel, sol] = potgetlpres(u, scal, 1);

pdratio = scal.bnrm / scal.cnrm;

if params.omega
    omega = rp / rd;
    omega = omega / pdratio;
else
    omega = 1.0;
end % End if

% Initialize parameters
beta = 0.1;
ncone = size(M, 2) - m + length(ubsub);
rho = params.rho * (ncone + sqrt(ncone));
potold = inf;

isok = potlogger(0, beta, omega, rp, rd, sol.pobj, sol.dobj,...
                  presrel, dresrel, gaprel, params.tol);

% Intialize projection
C = [z']; %#ok
CCT = C * C';
CCTinv = full(inv(CCT));

% Initialize iterations
uold = u;
z = u;
failed = false;

% Ready to go for potential reduction
for iter = 1:params.maxiter
    
    if mod(iter, 100) == 1 && params.omega
        potold = inf;
        omega = 0.4 * omega + 0.6 * (rp / rd) / pdratio;
%         beta = min(0.5, params.betamax);
    end % End if 
    
    if isok
        break;
    end % End if
    
    % Compute gradient
    [f, rp, rd, g] = potresi(m, M, u, v, omega);
    
    slu = ubsub - u(ubidinall);
    
    % Prepare gradient and momentum
    grad = g;
    grad(coneid) = grad(coneid) - (f ./ u(coneid)) / rho;
    grad(ubidinall) = grad(ubidinall) + (f./slu) / rho;
    mmtm = u - uold;

    % Gradient projection
    grad = grad - C' * (CCTinv * (C * grad));
    mmtm = mmtm - C' * (CCTinv * (C * mmtm));
    
    % Prepare hessian
    nrmgrad = norm(grad);
    nrmmmtm = norm(mmtm);
    
    if nrmmmtm > 0.0
        mmtm = mmtm / nrmmmtm;
    end % End if
    
    grad = grad / nrmgrad;
    nrmgrad = nrmgrad * rho / f;
    Mgrad = M * grad;
    Mmmtm = M * mmtm;
    
    gTgrad = g' * grad;
    gTmmtm = g' * mmtm;
    z(coneid) = 1./(u(coneid).^2);
    z(ubidinall) = z(ubidinall) + 1./(slu.^2);
    z(coneid) = sqrt(z(coneid));
    zinvgrad = grad .* z;
    zinvmmtm = mmtm .* z;
    zinvgrad(1:m) = 0;
    zinvmmtm(1:m) = 0;
    
    gradZZgrad = norm(zinvgrad)^2;
    mmtmZZmmtm = norm(zinvmmtm)^2;
    gradZZmmtm = zinvgrad' * zinvmmtm;
    
    gradHgrad = -rho * (gTgrad / f)^2 + gradZZgrad +...
                    rho * quadformw(m, omega, Mgrad, Mgrad) / f;
    mmtmHmmtm = -rho * (gTmmtm / f)^2 + mmtmZZmmtm +...
                    rho * quadformw(m, omega, Mmmtm, Mmmtm) / f;
    mmtmHgrad = -rho * (gTgrad * gTmmtm) / f^2 + gradZZmmtm +...
                    rho * quadformw(m, omega, Mmmtm, Mgrad) / f;
    
    gradTgrad = grad' * grad * nrmgrad;
    gradTmmtm = grad' * mmtm * nrmgrad;
    
    H = [gradHgrad, mmtmHgrad;
         mmtmHgrad, mmtmHmmtm];
    h = [gradTgrad; gradTmmtm];
    
    G = [gradZZgrad, gradZZmmtm;
         gradZZmmtm, mmtmZZmmtm];
     
    % Solve trust region subproblem
    while true
        
        [alpha, potrdcm] = subtrust(H, h, G, beta^2 / 2.5, 1e-10);
        d = alpha(1) * grad + alpha(2) * mmtm;
        d = d - C' * (CCTinv * (C * d));
        
        % Measure true potential reduction
        utmp = u + d;
        
        if (min(utmp(coneid)) < 0) || (min([ubsub - utmp(ubidinall); 0]) < 0)
            beta = beta * 0.5;
            continue;
        end % End if
        
        [ftmp, ~, ~, ~] = potresi(m, M, utmp, v, omega);
        potnew = rho * log(ftmp) - sum(log(utmp(coneid))) -...
                    sum(log(ubsub - utmp(ubidinall)));
                
        potrdc = potnew - potold;
        ratio = potrdc / potrdcm;
        
        if beta < 1e-05
            failed = true;
            break;
        end % End if
        
        % Adaptive trust region size
        if (ratio < 0.25 || potrdc > 0.0)
            beta = beta * 0.5;
        else
            if ratio > 0.75
                beta = min(beta * 2, params.betamax);
            elseif ratio < 0.5
                beta = beta * 0.8;
            end % End if
            uold = u;
            u = u + d;
            potold = potnew;
            break;
        end % End if
    end % End while
    
    if mod(iter, params.log) == 1 || (params.log == 1)
        % Evaluate new point
        [~, ~, presrel, dresrel, ~, gaprel, sol] = potgetlpres(u, scal, 1);
        isok = potlogger(iter, beta, omega, rp, rd, sol.pobj, sol.dobj,...
                        presrel, dresrel, gaprel, params.tol);
    end % End if
    
    if failed
        fprintf("Algorithm fails due to invalid trust radius \n");
        break;
    end % End if
end % End for

x = sol.x;
y = sol.y;
s = sol.s;
w = sol.w;

end % End function