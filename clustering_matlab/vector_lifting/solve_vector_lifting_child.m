function [result] = solve_vector_lifting_child(data, card, ML, CL, init_B_cell, global_ub, global_X, params)

    % SOLVE VECTOR-LIFTING SDP FOR ccMSSC no shrinking

    % data: data matrix (n x d)
    % card: cardinality matrix (k x k)
    % ML: matrix of must-link constraints (n_ml x 2)
    % CL: matrix of cannot-link constraints (n_cl x 2)
    % init_B_cell: cell array of initial inequalities (for empty cell pass cell(1, k))
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
    % params.cp_tol: cutting-plane tolerance (default 1e-4)
    % params.cp_percineq: percentage of inequalities to add (0.1 = 10%)
    % params.cp_epsineq: tolerance for checking the violation of inequalities (default 1e-4)
    % params.cp_activeineq: tolerance for checking active inequalities (default 1e-5) 
    % params.cp_inheritineq: inherit inequalitis from paret node (default yes - 1)
    % params.gurobi_vebose: verbosity level of gurobi (default 0)
    
    warning off;
    
    %disp(params);
    
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
    [blk, At, C, b, L, B_cell_cl, ~, n_cl_block] = build_vector_lifting_pairwise(data, card, ML, CL);
    
    init_n_ineq_block = zeros(k, 1);
    for j=1:k
        init_n_ineq_block(j) = size(init_B_cell{j}, 2);
    end
    Bt = cell(1, k);
    l = [];
    for j=1:k
        B = [Y_zeros_cell(sum(init_n_ineq_block(1:(j-1)))+sum(n_cl_block(1:(j-1))), n+1), ... 
            init_B_cell{j}, B_cell_cl{j}, ... 
            Y_zeros_cell(sum(init_n_ineq_block(j+1:k))+sum(n_cl_block(j+1:k)), n+1)];
        l = [l; zeros(init_n_ineq_block(j), 1); -ones(n_cl_block(j), 1)];
        Bt(j) = svec(blk(j, :), B, 1);
    end
    
    options.printlevel = params.sdp_verbose;
    options.tol = params.sdp_tol;
    options.stopoption = 0;
    options.AATsolve.method = 'iterative';
    options.maxiter = params.sdp_maxiter;
    options.maxtime = params.sdp_maxtime;

    % solve SDP
    [~, Yopt, ~, y, Z1, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
    
    % recover xj blocks
    X = zeros(n, k);
    PI_cell = cell(1, k);
    for j=1:k
        Yj = Yopt{j}; % get j-th block
        X(:, j) = Yj(2:n+1, 1);
        PI_cell{j} = Yj(2:n+1, 2:n+1);
    end

    % post-processing
    [lb0, lb1] = safe_bound_error_vl(blk, At, C, b, y, Z2, Bt, y2, l, Yopt, card);
    [~, lb2] = safe_bound_lp_vl(blk, Z1, At, Bt, C, b, l, params.n_threads);
    % safe lower bound
    lb = max(lb1, lb2);
    
    % heuristic with the SDP solution
    [ub_X, Xass_X] = cardinality_heuristic_pairwise(X, [], data, card, ML, CL, params.gurobi_verbose);
    [ub_Z, Xass_Z] = cardinality_heuristic_pairwise([], PI_cell, data, card, ML, CL, params.gurobi_verbose);
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
    result.ineq_list = [result.ineq_list; sum(init_n_ineq_block)];
    
    result.best_lb = lb;
    if ub <= global_ub
        result.best_ub = ub;
        result.best_Xass = sparse(Xass);
    else
        result.best_ub = global_ub;
        result.best_Xass = sparse(global_X);
    end
    result.best_Xlb = X; % best SDP feasible assignment [0, 1] matrix
    result.best_PI = PI_cell; % best SDP feasibile PI k blocks nxn
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
    l = [];
    n_ineq_block = init_n_ineq_block;

    for i=1:params.cp_maxiter
        
        if ~isempty(l)
           
            fprintf('Inequalities before removal = %d: ', sum(n_ineq_block));
            for j=1:k
                fprintf('%d ', n_ineq_block(j));
                Yvec = svec(blk(j, :), Yopt{j}, 1); 
                Bvec = svec(blk(j, :), B_cell{j}, 1);
                active_idx = abs(Bvec{1}' * Yvec) <= params.cp_activeineq;
                temp_B = B_cell{j};
                B_cell{j} = temp_B(active_idx);
                n_active_j = size(B_cell{j}, 2);
                n_ineq_block(j) = n_active_j;
            end
            %l = zeros(sum(n_ineq_block), 1);
            fprintf('\nInequalities after removal = %d: ', sum(n_ineq_block));
            for j=1:k
                fprintf('%d ', n_ineq_block(j));
            end
            fprintf('\n');
            
        end
                
        % separate valid inequalities
        separated_ineq = 0;
        for j=1:k
            [temp_Bcell, n_ineq] = separate_triangle_P(Yopt{j}, ...
                params.cp_epsineq, params.cp_maxineq, params.cp_percineq*(1/k));
            B_cell{j} = [B_cell{j}, temp_Bcell];
            n_ineq_block(j) = n_ineq_block(j) + n_ineq;
            %l = [l; zeros(n_ineq, 1)];
            separated_ineq = separated_ineq + n_ineq;
        end
        
        if separated_ineq <= 1
            result.cp_flag = 0;
            break
        end
        
        Bt = cell(1, k);
        l = [];
        for j=1:k
            B = [Y_zeros_cell(sum(n_ineq_block(1:(j-1)))+sum(n_cl_block(1:(j-1))), n+1), ...
                B_cell{j}, B_cell_cl{j}, ...
                Y_zeros_cell(sum(n_ineq_block(j+1:k))+sum(n_cl_block(j+1:k)), n+1)];
            l = [l; zeros(n_ineq_block(j), 1); -ones(n_cl_block(j), 1)];
            Bt(j) = svec(blk(j, :), B, 1);
        end
                
        % print inequalities added for each block j=1,...,k
        tot_ineq = sum(n_ineq_block);
        fprintf('Inequalities added = %d: ', tot_ineq);
        for j=1:k
            fprintf('%d ', n_ineq_block(j));
        end
        fprintf('\n');
        
        % solve SDP with inequalities
        [~, Yopt, ~, y, Z1, Z2, y2, ~, info, ~] = sdpnalplus(blk, At, C, b, L, [], Bt, l, [], options);
        result.cp_iter = result.cp_iter + 1;
        
        % SDP is not solved successfully
        if info.termcode == 1
            result.cp_flag = -1;
            break
        end
        
        % recover xj blocks
        X = zeros(n, k);
        PI_cell = cell(1, k);
        for j=1:k
            Yj = Yopt{j};
            X(:, j) = Yj(2:n+1, 1);
            PI_cell{j} = Yj(2:n+1, 2:n+1);
        end
        % post-processing
        [lb0, lb1] = safe_bound_error_vl(blk, At, C, b, y, Z2, Bt, y2, l, Yopt, card);
        [~, lb2] = safe_bound_lp_vl(blk, Z1, At, Bt, C, b, l, params.n_threads);
        lb = max(lb1, lb2);
        
        % heuristic with the SDP solution
        [ub_X, Xass_X] = cardinality_heuristic_pairwise(X, [], data, card, ML, CL, params.gurobi_verbose);
        [ub_Z, Xass_Z] = cardinality_heuristic_pairwise([], PI_cell, data, card, ML, CL, params.gurobi_verbose);
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
        result.ineq_list = [result.ineq_list; tot_ineq];
        
        if result.best_ub - ub >= 1e-6
            % update best upper bound
            result.best_ub = ub;
            result.best_Xass = sparse(Xass);
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
            result.best_PI = PI_cell;
            gap = (result.best_ub - result.best_lb) / result.best_ub;
            result.gap_list = [result.gap_list; gap];
            result.cp_flag = 3;
            result.best_B_cell = B_cell;
            if gap <= params.bb_tol
                result.cp_flag = 1;
            end
            break
        end
        
        % update best lower bound
        result.best_lb = lb;
        result.best_Xlb = X;
        result.best_PI = PI_cell;
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
        result.best_B_cell = cell(1, k);
    end

end
