function [alphabest] = adsqp(AX, Xc, target)
% Solve the QP from anderson acceleration
% minimize    0.5 * ||AXa||^2
% subject to   Xc * a >= 0
%              e' * a  = 1

Q = AX' * AX;
V = Xc;

[nc, k] = size(Xc);

alpha = ones(k, 1) / k;
nu = 0.0;
s = ones(nc, 1);
lbd = ones(nc, 1);
e = ones(k, 1);
sigma = 0.1;

pobjold = inf;

pbest = pobjold;
alphabest = alpha;

l = ones(nc, 1) * 0.0;
l(end) = min(Xc(end, :));

for i = 1:100
    
    mu = s' * lbd / nc;
    sinvl = lbd ./ s;
    SLV = diag(sqrt(sinvl)) * V;
    VTSLV = SLV' * SLV;
    sinv = s.^-1;
    M = Q + VTSLV + eye(k) * 1e-15;
    r = M * alpha + e * nu - V' * (lbd + mu * sigma * sinv) - V' * (sinvl .* l);
    
    M = sparse(M);
    d1 = cholmod2(M, e);
    d2 = cholmod2(M, r);
%     decomp = decomposition(M, 'chol');
%     d1 = decomp \ e;
%     d2 = decomp \ r;
    
    dnu = (sum(alpha) - 1 - e' * d2) / (e' * d1);
    dalpha = -d1 * dnu - d2;
    Vaplusda = V * (alpha + dalpha);
    dlbd = mu * sigma * sinv - sinvl .* Vaplusda;
    ds = Vaplusda - s - l;
    
    % Ratio test
    step = 0.995 / abs(min([dlbd./lbd; ds./s]));
    step = min(step, 1);
    s = s + step * ds;
    lbd = lbd + step * dlbd;
    nu = nu + step * dnu;
    alpha = alpha + step * dalpha;
    pobj = 0.5 * norm(AX * alpha)^2;
    
    if pobj < pbest
        alphabest = alpha;
    end % End if
    
    inf1 = V * alpha - s - l;
    inf2 = Q * alpha - V' * lbd + nu * e;
    fprintf("%3d %10.8e %10.3e %10.3e %10.3e %10.3e\n", i, pobj, abs(e' * alpha - 1), norm(inf1), norm(inf2), mu);
    
    if mu < 1e-07
        break;
    end % End if
    
%     if (pobjnew < pobjold) && (pobjold - pobjnew < 1e-06 * pobj)
%         break;
%     end % End if

end % End for

if isnan(pobj)
    keyboard;
end % End if

end % End function