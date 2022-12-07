function [] = test_netlib(fname, param, maxmn, minmn, fileID)

data = preprocess(fname);
A = data.A;
b = data.b;
c = data.c;

[m, n] = size(A);

if max(m, n) > maxmn || min(m, n) < minmn
    fprintf(fileID, "| %30s | %+3.1e | %+3.1e | %3.1e | %3.1e | %3.1e | %5.1f | %s \n",...
        fname, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, "Ignored");
    return;
end % End if

tic;

try
    [x, y, s] = potlp(A, b, c, param);
catch
    x = zeros(n, 1);
    s = zeros(n, 1);
    y = zeros(m, 1);
end % End try

t = toc;
pobj = c' * x;
dobj = b' * y;
pres = norm(A * x - b) / ( 1 + norm(b, 1) );
dres = norm(A' * y + s - c) / ( 1 + norm(c, 1) );
cpl = abs(pobj - dobj) / (1 + abs(pobj) + abs(dobj));
max1 = max([pres, dres, cpl]);
mm = max1;

if mm < 1e-04
    status = "Optimal";
elseif mm < 1e-03
    status = "Inaccurate";
else
    status = "Failed";
end % End if

fprintf(fileID, "| %30s | %+3.1e | %+3.1e | %3.1e | %3.1e | %3.1e | %5.1f | %s \n",...
    fname, pobj, dobj, pres, dres, cpl, t, status);

end % End function