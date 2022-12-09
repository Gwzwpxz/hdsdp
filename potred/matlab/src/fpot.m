function [fval, grad, resi] = fpot(A, ATA, x)

grad = ATA * x;
resi = A * x;
nrm = norm(resi);
fval = 0.5 * nrm^2;

end % End function