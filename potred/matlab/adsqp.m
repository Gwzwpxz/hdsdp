function [alpha] = adsqp(AX, Xc)
% Solve the QP from anderson acceleration
% minimize    0.5 * ||AXa||^2
% subject to   Xc * a >= 0
%              e' * a  = 1

Q = AX' * AX;
V = Xc;

[nc, k] = size(Xc);

alpha = zeros(k, 1);
alpha(end) = 1.0;
nu = 0.0;
s = ones(nc, 1);
lbd = ones(nc, 1);
e = ones(k, 1);
sigma = 0.8;

pobjold = inf;
pobjnew = pobjold;

for i = 1:100
    
    mu = s' * lbd / nc;
    VTSLV = V' * diag(lbd ./ s) * V;
    sinv = s.^-1;
    r = Q * alpha - V' * lbd + e * nu - mu * sigma * V' * sinv + VTSLV * alpha;
    M = Q + VTSLV;
    
    d1 = M \ e;
    d2 = M \ r;
    dnu = -(e' * d2) / (e' * d1);
    dalpha = -d1 * dnu - d2;
    dlbd = mu * sigma * sinv - ((V * (alpha + dalpha)) ./ s) .* lbd;
    ds = V * dalpha + (V * alpha - s);
    
    % Ratio test
    step = 0.99995 / abs(min([dlbd./lbd; ds./s]));
    step = min(step, 1);
    s = s + step * ds;
    lbd = lbd + step * dlbd;
    nu = nu + step * dnu;
    alpha = alpha + step * dalpha;
    pobj = 0.5 * norm(AX * alpha)^2;
    
    inf1 = V * alpha - s;
    inf2 = Q * alpha - V' * lbd + nu * e;
    fprintf("%3d %10.6e %10.3e %10.3e %10.3e %10.3e\n", i, pobj, abs(e' * alpha - 1), norm(inf1), norm(inf2), mu);
    
    pobjold = pobjnew;
    pobjnew = pobj;
    
    if (pobjnew < pobjold) && (pobjold - pobjnew < 1e-06 * pobj)
        break;
    end % End if

end % End for

end % End function