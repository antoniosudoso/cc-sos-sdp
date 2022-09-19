function [ub, Xopt] = cardinality_heuristic_pairwise(Xbar, Zbar, Ws, C, ML, CL, verbose)

    % Xbar: solution of the SDP relaxation n x k
    % Zbar: solution of the SDP relaxation n x n
    % Ws: data matrix n x d 
    % C: cardinality matrix k x k
    % ML: must-link constraints matrix n_ml x 2
    % CL: cannot-link constraints matrix n_cl x 2
    % verbose: display or not gurobi log
    
    [n, d] = size(Ws);
    k = size(C, 1);
    
    n_ml = size(ML, 1);
    n_cl = size(CL, 1);
    
    % Index helper function
    assidx = @(i, j) i+(j-1)*n;
    
    % Build model
    model.vtype = repmat('B', n*k, 1);
    model.A = sparse(n+k+k*n_ml+k*n_cl, n*k);
    count = 1;
    for i=1:n
        for j=1:k
            model.A(i, assidx(i, j)) = 1;
        end
        count = count + 1;
    end
    for j=1:k
        model.A(j+n, 1+n*(j-1):n*j) = 1;
        count = count + 1;
    end
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        for h=1:k
            model.A(count, assidx(i, h)) = 1;
            model.A(count, assidx(j, h)) = -1;
            count = count + 1;
        end
    end
    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        for h=1:k
            model.A(count, assidx(i, h)) = 1;
            model.A(count, assidx(j, h)) = 1;
            count = count + 1;
        end
    end
    % spy(model.A)
    model.rhs = [ones(n, 1); diag(C); zeros(k*n_ml, 1); ones(k*n_cl, 1)]; % n + k + k*n_ml + k*n_cl constraints
    model.sense = [repmat('=', n+k+k*n_ml, 1); repmat('<', k*n_cl, 1)];
        
    
    if isempty(Xbar) && isempty(Zbar)
        % random cluster centers
        rng(1727);
        index = randsample(1:n, k);
        init_center = Ws(index, :);
    end
    
    if ~isempty(Xbar) 
        %fprintf('Finding nearest assignment matrix...\n');
        model.modelsense = 'max';
        model.obj = Xbar(:);
        params.outputflag = verbose;
        % params.threads = n_threads;
        result = gurobi(model, params);
        % disp(result.objval);
        init_X = reshape(result.x, [n, k]);
        init_center = zeros(k, d);
        for j=1:k
            [~, q] = max(init_X, [], 2);
            init_center(j, :) = sum(Ws(q == j, :)) / C(j, j);
        end
    end
    
    if ~isempty(Zbar)
        new_Zbar = zeros(n, n);  
        if iscell(Zbar)
            for j=1:k
                new_Zbar = new_Zbar + (1/C(j,j)) * Zbar{j};
            end 
        else
            new_Zbar = Zbar;
        end
        %fprintf('Finding best rank-k approximation...\n');
        % spectral decomposion of psd matrix
        [Ubar, Dbar] = eig(new_Zbar);
        [~, ind] = sort(diag(Dbar));
	    Dbar = Dbar(ind, ind);
	    Ubar = Ubar(:, ind);
        % rank k approximation by truncating the spectral decomposition
        new_Zopt = Ubar(:, n-k+1:n)*Dbar(n-k+1:n, n-k+1:n)*Ubar(:, n-k+1:n)';
        M = new_Zopt * Ws;
        rng(1727);
        [~, init_center, ~, ~] = kmeans(M, k, 'Start', 'plus', 'Replicates', 100);
    end
    
    ub = inf;
    Xopt = [];
    max_iter = 100;
    trace_const = trace(Ws * Ws');
    params.outputflag = verbose;
    % params.threads = n_threads;
    for h=1:max_iter
        % Build model
        model.modelsense = 'min';
        obj = Ws * init_center';
        model.obj = -2 * obj(:);
        result = gurobi(model, params);
        current_Xopt = reshape(result.x, [n, k]);
        obj_const = trace(init_center * init_center' * C);
        current_ub = result.objval + trace_const + obj_const;
        if abs(current_ub - ub) <= 1e-7
            break
        else
            % update upper bound and solution
            ub = current_ub;
            Xopt = current_Xopt;
        end
        %fprintf('Found upper bound = %10.9e \n', ub);
        % compute centroids
        init_center = zeros(k, d);
        for j=1:k
            [~, q] = max(current_Xopt, [], 2);
            init_center(j, :) = sum(Ws(q == j, :)) / C(j, j);
        end
    end
     
end
