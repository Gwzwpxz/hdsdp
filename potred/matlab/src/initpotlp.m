function [u, z] = initpotlp(m, n, M, v, ubidinall, ubsub, coneid, scal) %#ok
% Initialize the potential solver
% Upper-unbounded conic variables are initialized at 1.0
% Other conic variables take ub / 2
% Non-conic variables are initialized to ensure that the duality gap is 0

v; %#ok
u = gdinit(M, v, m, 1, 200);
u(coneid) = 1.0; % max(u(coneid), 1000);
% u(m+1:m+n) = norm(v(1:m));
% u(m+n+1:end) = norm(v(m+1:end));
u(ubidinall) = min(ubsub / 2);

if scal.method == 'A'
    % Interior point warm start
%     A = scal.A;
%     b = scal.b;
%     c = scal.c;
%     [x, ~, s, ~, tau] = hsdipm(A, b, c);
%     u(scal.xid) = scal.E .* (x / tau);
%     u(scal.sid) = (s / tau) ./ scal.E;
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