function [isok] = potlogger(iter, beta, omega, rp, rd, pobj, dobj, pres, dres, gap, tol) %#ok
% Logging for potential reduction

if iter == 0
    fprintf("%4s %8s %8s %8s %8s %8s [%5s, %4s]\n",...
        "iter", "pobj", "dobj", "gap", "pres", "dres", "omega", "beta");
end % End if

isok = false;
fprintf("%4d %+8.1e %+8.1e %8.1e %8.1e %8.1e [%3.3f, %2.1f]\n", iter, pobj, dobj, gap, pres, dres, omega, beta);
if pres < tol && dres < tol && gap < tol
    isok = true;
end % End if

end % End function