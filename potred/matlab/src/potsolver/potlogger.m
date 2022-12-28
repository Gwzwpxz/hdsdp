function [isok] = potlogger(iter, omega, rp, rd, pobj, dobj, pres, dres, gap, tol) %#ok
% Logging for potential reduction

if iter == 0
    fprintf("%4s %8s %8s %8s %8s %8s [%3s]\n",...
        "iter", "pobj", "dobj", "gap", "pres", "dres", "omg");
end % End if

isok = false;
fprintf("%4d %+8.1e %+8.1e %8.1e %8.1e %8.1e [%3.3f]\n", iter, pobj, dobj, gap, pres, dres, omega);
if pres < tol && dres < tol && gap < tol
    isok = true;
end % End if

end % End function