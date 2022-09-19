function [y, LB] = safe_bound_lp_ml(blk, Z, At, Bt, C, b, l, threads)

    % Determine a feasible y, given Z and S
    % Dual SDP: max b'y s.t. Aty + Z + S = C;
    % Z psd, S nonnegative
    
    %fprintf('\nRunning LP post-processing... \n');
    
    n = size(C, 1);
    n2 = n*n;
    if n2 == size(At, 2)
        At = At';
    end

    if n2 == size(Bt, 2)
        Bt = Bt';
    end

    [V, L] = eig(Z); 
    idx = find(L < 0); 
    numneg = length(idx); 
    if numneg
        % project Z onto psd cone
        Zplus  = V * max(zeros(n, n), L) * V';
    else 
        % Z is psd
        Zplus = Z;
    end

    %rhs
    M = C - Zplus;
    M = svec(blk, M, 1);

    options.Display = 'off';
    options.Method = 'barrier'; % dual simplex, barrier, concurrent
    options.Threads = double(threads);
    if ~isempty(Bt) && ~isempty(l)
        % with inequalities
        [y, LB1, flag] = my_linprog([-b; -l], [At, Bt], M, [], [], ...
            [-inf * ones(length(b), 1); zeros(length(l), 1)], ...
            inf * ones(length([b; l]), 1), options);
    else
        % without inequalities
        [y, LB1, flag] = my_linprog(-b, At, M, [], [], [], [], options);
    end

    if flag ~= 1
        LB = inf;
    else
        LB = -LB1;
    end

    %fprintf('LP safe lower bound = %10.9e \n', LB);
end
