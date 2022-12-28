function [scal, M, v, ubid, ubsub, coneid] = ...
    lpscale(A, b, c, ub, ruiz, l2scal, method)
%{
Apply
    a. L-inf  Ruiz scaling
    b. L-2    Row/column matrix equilibration
    c. p-norm Scaling of b and c
to improve conditioning of LP.

There are two ways to do the scaling
Method M.
    Scaling is performed after setting up the M matrix
On exit,
    a. M = D [ 0  A  0   0 ] E
             [ A' 0  I -I_w]
    b. v = D [b; c]
    c. ubs = ub .* e_x
    
Method A
    Scaling is performed on A matrix and M is set up with scaled A
On exit,
    a. A <- DAE
    b. M =  [ 0  A  0   0 ]
            [ A' 0  I -I_w]
    c. v = [ D * b; E * c]
    c. ubs = ub .* e_x
    
    %}
    
    % Retrieve dimension
    [m, n] = size(A);
    
    if nargin < 7
        % Using Method A by default
        method = 'A';
    end % End
    
    % Common matrix components
    ubid = find(ub < +inf);
    nub = length(ubid);
    Iw = sparse(ubid, 1:nub, ones(nub, 1), n, nub);
    
    scal.method = method;
    scal.A = A;
    scal.b = b;
    scal.c = c;
    scal.ub = ub;
    scal.ubid = ubid;
    scal.bnrm = norm(b);
    scal.cnrm = norm(c);

    yid = 1:m;
    xid = m+1:m+n;
    sid = m+n+1:m+n+n;
    wid = m+n+n+1:m+n+n+length(ubid);
    
    scal.yid = yid;
    scal.xid = xid;
    scal.sid = sid;
    scal.wid = wid;
    
    fprintf("LP statistics: ||b|| = %5.3e. ||c|| = %5.3e \n", scal.bnrm, scal.cnrm);
    fprintf("Using scaling Method %s \n", method);
    
    if method == 'M'
        % Build M and scale
        M = [sparse(m, m), A,            sparse(m, n + nub);
             A',           sparse(n, n), speye(n), -Iw];
        
        % Scale matrix
        [Druiz, Eruiz, M] = ruizscale(M, ruiz);
        if l2scal
            [Dl2, El2, M] = l2scale(M);
            D = full(Druiz .* Dl2);
            E = full(Eruiz .* El2);
        else
            D = full(Druiz);
            E = full(Eruiz);
        end % End if
        
        % Scale bounds
        ubsub = ub(ubid) .* E(m+ubid);
        
        % Scale primal and dual objective
        v = [b; c];
        v = v ./ D;
    else
        % Scale A and build M
        [Druiz, Eruiz, A] = ruizscale(A, ruiz);
        if l2scal
            [Dl2, El2, A] = l2scale(A);
            D = full(Druiz .* Dl2);
            E = full(Eruiz .* El2);
        else
            D = full(Druiz);
            E = full(Eruiz);
        end % End if
        
        % Scale bounds
        ubsub = ub(ubid) .* E(ubid);
        
        % Scale primal and dual objective
        v = [b./D; c./E];
        
        % Build M
        M = [sparse(m, m), A,            sparse(m, n + nub);
             A',           sparse(n, n), speye(n), -Iw];
        
    end % End if
    
    % Save scaling
    scal.D = D;
    scal.E = E;
    
    % Scatter ub
    coneid = (m + 1):size(M, 2);
    ubid = coneid(ubid);
    
end % End function