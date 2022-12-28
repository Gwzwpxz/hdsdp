function [Px] = proj2d(n, sumc, nrmb, nrmc, C, x)
% Projection onto joint Null space
CCTinv = [nrmb^2 + nrmc^2, sumc;
          sumc, 2 * n + 1];
CCTinv = CCTinv / ((2 * n + 1) * (nrmb^2 + nrmc^2) - sumc^2);
Px = x - C' * (CCTinv * (C * x));
end % End function