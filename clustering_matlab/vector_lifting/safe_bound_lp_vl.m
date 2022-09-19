function [y, LB] = safe_bound_lp_vl(blk, Z, At, Bt, C, b, l, threads)

    % Determine a feasible y, given Z and S
    % Dual SDP: max b'y s.t. Aty + Z + S = C;
    % Z psd, S nonnegative
    
    %fprintf('\nRunning LP post-processing... \n');

    n = blk{1, 2};
    k = length(blk);
    At_list = [];
    Bt_list = [];
    rhs_list = [];

    for j=1:k

        [V, L] = eig(Z{j}); 
        idx = find(L < 0); 
        numneg = length(idx); 
        if numneg
            % project Z onto psd cone
            Zplus  = V * max(zeros(n, n), L) * V';
        else 
            % Z is psd
            Zplus = Z{j};
        end
        %rhs
        M = C{j} - Zplus;
        M = svec(blk(j, :), M, 1);
        rhs_list = [rhs_list; M];
        At_list = [At_list; At{j}];
        if ~isempty(Bt)
            Bt_list = [Bt_list; Bt{j}];
        end

    end
    
    options.Display = 'off';
    options.Method = 'barrier';
    options.Threads = double(threads);
    if ~isempty(Bt) && ~isempty(l)
        % with inequalities
        [y, LB1, flag] = my_linprog([-b; -l], [At_list, Bt_list], rhs_list, [], [], ...
            [-inf * ones(length(b), 1); zeros(length(l), 1)], ...
            inf * ones(length([b; l]), 1), options);
    else
        % without inequalities
        [y, LB1, flag] = my_linprog(-b, At_list, rhs_list, [], [], [], [], options);
    end

    if flag ~= 1
        LB = inf;
    else
        LB = -LB1;
    end

    %fprintf('LP safe lower bound = %10.9e \n', LB);
end
