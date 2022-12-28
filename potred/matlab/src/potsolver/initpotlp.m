function [u, z] = initpotlp(m, n, M, v, ubidinall, ubsub, coneid, scal)
% Initialize the potential solver
% Upper-unbounded conic variables are initialized at 1.0
% Other conic variables take ub / 2
% Non-conic variables are initialized to ensure that the duality gap is 0

v; %#ok
u = zeros(size(M, 2), 1);

u(coneid) = 1;
u(ubidinall) = ubsub / 2;

if scal.method == 'A'
    vp = v(1:m);
    vd = v(m+1:end);
    alpha = vd' * u(scal.xid) + ubsub' * u(scal.wid);
    u(scal.yid) = alpha * vp / norm(vp)^2;
    z = [-vp; vd; zeros(n, 1); ubsub];
else
    utmp = u ./ scal.E; 
    alpha = scal.c' * utmp(scal.xid) + scal.ub(scal.ubid)' * utmp(scal.wid);
    utmp(scal.yid) = alpha * scal.b / norm(scal.b)^2;
    u = utmp .* scal.E;
    z = [-scal.b; scal.c; zeros(n, 1); scal.ub(scal.ubid)];
    z = z ./ scal.E;
end % End if

end % End function