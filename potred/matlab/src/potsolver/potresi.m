function [fval, rp, rd, grad] = potresi(m, M, u, v, omega)
% Compute gradient of re-parametrized objective
% f_w(u) = 0.5 ||M * u - v||^2_w

r = M * u - v;
rp = norm(r(1:m));
rd = norm(r(m+1:end));
fval = omega * rp^2 + (1/omega) * rd^2;
fval = fval / 2;

if omega ~= 1.0
    r(1:m) = r(1:m) * omega;
    r(m+1:end) = r(m+1:end) / omega;
end % End if
grad = M' * r;

end % End function