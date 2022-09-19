function [result] = solve_matrix_lifting_child(data, card, ML, CL, init_B_cell, global_ub, global_X, params)

    % SOLVE MATRIX-LIFTING SDP FOR ccMSSC no shrinking
    
    % data: data matrix (n x d)
    % card: cardinality matrix (k x k)
    % ML: matrix of must-link constraints (n_ml x 2)
    % CL: matrix of cannot-link constraints (n_cl x 2)
    % init_B_cell: cell array of initial inequalities (for empty pass cell(0))
    % global_ub: best known upper bound (double)
    % global_X: best known assignment 0-1 matrix (n x k)
    
    % params.n_threads: number of threads for the current session
    % params.bb_tol: branch-and-bound tolerance (default 1e-3)
    % params.sdp_verbose: verbosity level of sdpnal+ (default 0)
    % params.sdp_tol: accuracy tolerance of sdpnal+ (default 1e-4)
    % params.sdp_maxiter: maximum number of iteration of sdpnal+ (default 50000)
    % params.sdp_maxtime: maximum solving time of sdpnal+ (default 3600)
    % params.cp_maxineq: maximum number of valid inequalities to separate (default 100000)
    % params.cp_maxiter: number of cutting-plane iterations (default 10)
    % params.cp_tol: cutting-plane tolerance (default 1e-5)
    % params.cp_percineq: percentage of inequalities to add (0.1 = 10%)
    % params.cp_epsineq: tolerance for checking the violation of inequalities (default 1e-4)
    % params.cp_activeineq: tolerance for checking active inequalities (default 1e-6)
    % params.cp_inheritineq: inherit inequalitis from paret node (default yes - 1)
    % params.gurobi_vebose: verbosity level of gurobi (default 0)
    
    warning off;
    
    maxNumCompThreads(params.n_threads);
 
    result = struct();
    result.lb0_list = [];
    result.lb_list = [];
    result.ub_list = [];
    result.gap_list = [];
    result.termcode_list = [];
    result.iter_list = [];
    result.time_list = [];
    result.ineq_list = [];

    n = size(data, 1);
    k = size(card, 1);
    [blk, At, C, b, L, B_cell_cl, l_cl] = build_matrix_lifting_pairwise(data, card, ML, CL);
    
    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.AATsolve.method = 'iterative';
    options.maxiter = params.sdp_maxiter;
    options.maxtime = params.sdp_maxtime;
    
    Bt = svec(blk, [init_B_cell, B_cell_cl], 1);
    init_l = zeros(size(init_B_cell, 2), 1);
    l = [init_l; l_cl];

    % solve SDP
    [~, Yopt, ~, y, Z1, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);

    Y = Yopt{1}; % get Y block
    X = Y(1:k, k+1:n+k)'; % get X block
    Z = Y(k+1:n+k, k+1:n+k); % get Z block
    % post-processing
    max_card = max(card(:));
    [lb0, lb1] = safe_bound_error_ml(blk, At, C, b, y, Z2, Bt, y2, l, Yopt, max_card);
    [~, lb2] = safe_bound_lp_ml(blk, Z1{1}, At{1}, Bt{1}, C{1}, b, l, params.n_threads);
    % safe lower bound
    lb = max(lb1, lb2);
    
    % heuristic with the SDP solution
    [ub_X, Xass_X] = cardinality_heuristic_pairwise(X, [], data, card, ML, CL, params.gurobi_verbose);
    [ub_Z, Xass_Z] = cardinality_heuristic_pairwise([], Z, data, card, ML, CL, params.gurobi_verbose);
    if ub_X <= ub_Z
        ub = ub_X;
        Xass = Xass_X;
    else
        ub = ub_Z;
        Xass = Xass_Z;
    end
    
    result.lb0_list = [result.lb0_list; lb0];
    result.lb_list = [result.lb_list; lb];
    result.ub_list = [result.ub_list; ub];
    result.termcode_list = [result.termcode_list; info.termcode];
    result.iter_list = [result.iter_list; info.iter];
    result.time_list = [result.time_list; info.totaltime];
    result.ineq_list = [result.ineq_list; size(init_B_cell, 2)];
    
    result.best_lb = lb;
    if ub <= global_ub
        result.best_ub = ub;
        result.best_Xass = sparse(Xass);
    else
        result.best_ub = global_ub;
        result.best_Xass = sparse(global_X);
    end
    result.best_Xlb = X; % best SDP feasible assignment [0, 1] matrix
    result.best_Z = Z; % best SDP clustering matrix [0, 1]
    result.best_B_cell = init_B_cell;
    result.cp_iter = 0;
    result.cp_flag = -2;
    
    % -2 - maximum number of iterations
    % -1 - SDP not solved or partially solved successfully
    %  0 - no violated inequalities
    %  1 - node must be pruned
    %  2 - lower bound smaller than the previous one
    %  3 - lower bound not sufficiently large
    
    gap = (result.best_ub - result.best_lb) / result.best_ub;
    fprintf('\nRelative gap = %10.9e \n', gap);
    result.gap_list = [result.gap_list; gap];

    if gap <= params.bb_tol
        result.cp_flag = 1;
        return
    end
    
    B_cell = init_B_cell;
    l = init_l;
    
    for i=1:params.cp_maxiter
        
        % remove inactive inequalities
        if ~isempty(B_cell) && ~isempty(l)
            fprintf('Inequalities before removal = %d \n', size(B_cell, 2));
            Yvec = svec(blk, Y, 1);
            Bt = svec(blk, B_cell, 1);
            active_idx = abs(Bt{1}' * Yvec) <= params.cp_activeineq;
            B_cell = B_cell(active_idx);
            fprintf('Inequalities after removal = %d \n', size(B_cell, 2));
        end
        
        % separate valid inequalities
        [temp_B_cell_1, n_ineq_1] = separate_pair_Y(Y, k, ...,
            params.cp_epsineq, params.cp_maxineq, params.cp_percineq);
        [temp_B_cell_2, n_ineq_2] = separate_triangle_Y(Y, k, ...,
            params.cp_epsineq, params.cp_maxineq, params.cp_percineq);
        fprintf('Inequalities added = %d \n', n_ineq_1+n_ineq_2);
        
        if n_ineq_1+n_ineq_2 <= 1
            result.cp_flag = 0;
            break
        end
        
        temp_B_cell = [temp_B_cell_1, temp_B_cell_2];
        B_cell = [B_cell, temp_B_cell];
        current_ineq = size(B_cell, 2);
        Bt = svec(blk, [B_cell, B_cell_cl], 1);
        l = [zeros(current_ineq, 1); l_cl];
        % solve SDP with inequalities
        [~, Yopt, ~, y, Z1, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
        result.cp_iter = result.cp_iter + 1;
        
        % SDP is not solved successfully
        if info.termcode == 1
            result.cp_flag = -1;
            break
        end
        
        Y = Yopt{1};
        X = Y(1:k, k+1:n+k)';
        Z = Y(k+1:n+k, k+1:n+k);
        [lb0, lb1] = safe_bound_error_ml(blk, At, C, b, y, Z2, Bt, y2, l, Yopt, max_card);
        [~, lb2] = safe_bound_lp_ml(blk, Z1{1}, At{1}, Bt{1}, C{1}, b, l, params.n_threads);
        lb = max(lb1, lb2);
        
        % heuristic with the SDP solution
        [ub_X, Xass_X] = cardinality_heuristic_pairwise(X, [], data, card, ML, CL, params.gurobi_verbose);
        [ub_Z, Xass_Z] = cardinality_heuristic_pairwise([], Z, data, card, ML, CL, params.gurobi_verbose);
        if ub_X <= ub_Z
            ub = ub_X;
            Xass = Xass_X;
        else
            ub = ub_Z;
            Xass = Xass_Z;
        end
        
        result.lb0_list = [result.lb0_list; lb0];
        result.lb_list = [result.lb_list; lb];
        result.ub_list = [result.ub_list; ub];
        result.termcode_list = [result.termcode_list; info.termcode];
        result.iter_list = [result.iter_list; info.iter];
        result.time_list = [result.time_list; info.totaltime];
        result.ineq_list = [result.ineq_list; current_ineq];
        
        if result.best_ub - ub >= 1e-6
            % update best upper bound
            result.best_ub = ub;
            result.best_Xass = sparse(Xass);
            result.best_Z = Z;
        end
        
        % current bound is worst than previous lower bound
        if lb - result.best_lb <= 1e-6
            result.cp_flag = 2;
            gap = (result.best_ub - result.best_lb) / result.best_ub;
            result.gap_list = [result.gap_list; gap];
            if gap <= params.bb_tol
                result.cp_flag = 1;
            end
            break
        end
        
        cp_stop = (lb - result.best_lb) / lb;
        % current bound is not sufficiently larger than previous lower bound
        if cp_stop <= params.cp_tol
            result.best_lb = lb;
            result.best_Xlb = X;
            result.best_B_cell = B_cell;
            gap = (result.best_ub - result.best_lb) / result.best_ub;
            result.gap_list = [result.gap_list; gap];
            result.cp_flag = 3;
            if gap <= params.bb_tol
                result.cp_flag = 1;
            end
            break
        end
        
        % update best lower bound
        result.best_lb = lb;
        result.best_Xlb = X;
        result.best_Z = Z;
        result.best_B_cell = B_cell;
        
        gap = (result.best_ub - result.best_lb) / result.best_ub;
        fprintf('\nRelative gap = %10.9e \n', gap);
        result.gap_list = [result.gap_list; gap];
        % prune the node
        if gap <= params.bb_tol
            result.cp_flag = 1;
            break
        end
                
    end
    
    if ~params.cp_inheritineq
        result.best_B_cell = cell(0);
    end

end
