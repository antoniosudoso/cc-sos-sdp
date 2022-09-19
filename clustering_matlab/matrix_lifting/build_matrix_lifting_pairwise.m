function [blk, At, C, b, L, B_cell, l] = build_matrix_lifting_pairwise(data, card, ML, CL)

    % data: data matrix n x d
    % card: cardinality matrix k x k
    % ML: matrix of must-link constraints (n_ml x 2)
    % CL: matrix of cannot-link constraints (n_cl x 2)
    
    % return problem data sdpnal+ format
    
    n = size(data, 1);
    k = size(card, 1);
    
    D = 0.5 * squareform(pdist(data, 'squaredeuclidean'));

    W_full = zeros(n+k, n+k);
    W_full(k+1:n+k, k+1:n+k) = -D;
    
    Y_tl = Y_slice_C(n, k);
    [X_rowsum, X_colsum] = Y_slice_X(n, k);
    [Z_rowsum, Z_diag] = Y_slice_Z(n, k, card);
    
    [A_cell_X_ml, b_X_ml, A_cell_Z_ml, b_Z_ml] = add_must_link_Y(n, k, ML);
    [B_cell_X_cl, l_X_cl, A_cell_Z_cl, b_Z_cl] = add_cannot_link_Y(n, k, CL);
    
    % first data point assigned to first cluster
    % A_cell_first = cell(1);
    % A_cell_first{1} = sparse([1+k, 1], [1, k+1], [0.5, 0.5], n+k, n+k);
    % b_first = 1;

    Acell = [Y_tl, X_rowsum, X_colsum, Z_rowsum, Z_diag, ...
        A_cell_X_ml, A_cell_Z_ml, A_cell_Z_cl];
    b_card = [];
    for i = 1:k
        for j = i:k
            b_card = [b_card; card(i, j)];
        end
    end
    b = [b_card; ones(n, 1); diag(card); ones(n, 1); zeros(n, 1); ...
        b_X_ml; b_Z_ml; b_Z_cl];
    
    blk = cell(1);
    blk{1,1} = 's';
    blk{1,2} = n+k;

    At = svec(blk, Acell, 1);

    C = cell(1);
    C{1} = -sparse(W_full);

    L = cell(1);
    L{1} = 0; % we want Y >= 0
    
    B_cell = B_cell_X_cl;
    l = l_X_cl;

end