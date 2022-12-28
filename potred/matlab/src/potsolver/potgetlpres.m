function [pres, dres, presrel, dresrel, gap, gaprel, sol] = potgetlpres(u, scal, checkcone)
% Recover infeasibilities in terms of the original LP
% pres = ||A * x - b||_2
% relpres = pres / (1 + ||b||_2)
% dres = ||A' * y + s - c||_2
% reldres = dres / (1 + ||c||_2)
% gap = |c' * x - b' * y|
% relgap = gap / (1 + |c' * x| + |b' * y|)
% 
% Two different scaling matrices induce different recovering procedures
%
% Upperbound is always strictly respected due to algorithm design
% Dual infeasibility may be improved using orthogonal projection

A = scal.A;
b = scal.b;
c = scal.c;
m = size(A, 1);

xid = scal.xid;
yid = scal.yid;
sid = scal.sid;
wid = scal.wid;

ubid = scal.ubid;
ub = scal.ub;
D = scal.D;
E = scal.E;
nrmb = scal.bnrm;
nrmc = scal.cnrm;

if scal.method == 'M'
    u = u./E;
    x = u(xid);
    y = u(yid);
    s = u(sid);
    w = u(wid);
else
    x = u(xid) ./ E;
    y = u(yid) ./ D;
    s = u(sid) .* E;
    w = u(wid) .* E(ubid);
end % End if

if checkcone
    assert(min(u(m+1:end)) >= 0.0);
    assert(all(x <= ub + 1e-12));
end % End if

% Compute residuals
pres = A * x - b;
dres = A' * y + s - c;
dres(ubid) = dres(ubid) - w;
pobj = c' * x + scal.objcon;
dobj = b' * y - ub(ubid)' * w + scal.objcon;
gap = pobj - dobj;

pres = norm(pres);
dres = norm(dres);
gap = abs(gap);

presrel = pres / (1 + nrmb);
dresrel = dres / (1 + nrmc);
gaprel = gap / (1 + abs(pobj) + abs(dobj));

sol.pobj = pobj;
sol.dobj = dobj;
sol.x = x;
sol.y = y;
sol.s = s;
sol.w = w;

end % End function