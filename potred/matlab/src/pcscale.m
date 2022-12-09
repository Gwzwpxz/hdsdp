function [D, E, Ascal] = pcscale(A, maxiter)
% Implement the PC scaling algorithm for general matrix

if nargin == 1
    maxiter = 100;
end % End if

[m, n] = size(A);
D = ones(m, 1);
E = ones(n, 1);

for i = 1:maxiter
    
    dR = sqrt(sum(abs(A), 2));
    dC = sum(abs(A), 1).^(-0.5);
    R = diag(dR);
    C = diag(dC);
    A = R \  (A * C);
    D = D ./ dR;
    E = E .* dC';
    
end % End for

Ascal = A;

end % End function