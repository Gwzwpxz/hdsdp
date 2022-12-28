function [alpha, mval] = subtrust(H, h, M, delta, tol)
% Solve the trust region subproblem and measure the reduction in the model

% If only one direction is supplied
if H(2, 2) == 0.0
    M11 = M(1, 1); 
    
    if H(1, 1) > 0 && (h(1) / H(1, 1))^2 * M11 < 2 * delta
        alpha(1) = -h(1) / H(1, 1);
        alpha(2) = 0;
        mval = 0.5 * alpha(1)^2 * H(1, 1) + h(1) * alpha(1);
        return;
    end % End if
        
    alpha11 = sqrt(2 * delta / M11);
    mval = 0.5 * alpha11^2 * H(1, 1) + h(1) * alpha11;
    if mval > 0
        alpha11 = - alpha11;
        mval = 0.5 * alpha11^2 * H(1, 1) + h(1) * alpha11;
    end % End if
    alpha(1) = alpha11;
    alpha(2) = 0.0;
    return;
end % End if

% If H is positive definite
try 
    chol(H);
    alpha = -H \ h;
    if alpha' * M * alpha <= 2 * delta
        mval = 0.5 * alpha' * H * alpha + h' * alpha;
        return;
    end % End if
catch
    
end % End try

% Otherwise dive into the subspace
L = chol(M)';
LinvQLTinv = L \ ((L \ H)');
lamlb = - eigs(LinvQLTinv, 1, 'smallestreal');

% Bisection solver
if lamlb < 0
    lamlb = 0.0;
end % End if

lamub = max(lamlb, 1);

while true
    lamub = lamub * 2;
    HplM = H + lamub * M;
    alpha = - HplM \ h;
    nrmdiff = alpha' * M * alpha - 2 * delta;
    if nrmdiff <= 0
        lamlb = lamub / 2;
        break;
    end % End if
end % End while

diff = lamub - lamlb;
    
while diff > 1e-04
    
    lam = lamlb + diff / 2;
    HplM = H + lam * M;
    alpha = - HplM \ h;
    nrmdiff = alpha' * M * alpha - 2 * delta;
    if nrmdiff > tol
        lamlb = lam; 
    elseif nrmdiff < -tol
        lamub = lam; 
    else
        break;
    end % End if
    diff = lamub - lamlb;
        
end % End while
       
mval = 0.5 * alpha' * H * alpha + h' * alpha;
if mval > 0
    mval = -1e+10;
end % End if
        
end % End function