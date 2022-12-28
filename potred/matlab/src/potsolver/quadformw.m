function [xWy] = quadformw(m, omega, x, y)
% Evaluate quadratic form x' * W * y

xWy = omega * (x(1:m)' * y(1:m)) + (x(m+1:end)' * y(m+1:end)) / omega;

end % End function