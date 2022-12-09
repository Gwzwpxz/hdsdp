function [pres, dres, cpl] = getresi(x, res, m, n)

pres = res(1:m);
dres = res(m + 1 : m + n + 1);
cpl = res(m + n + 1);

pres = pres / x(end);
dres = dres / x(end);
cpl = cpl / x(end);

end % End function