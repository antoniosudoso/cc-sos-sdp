function [blk, At, C, b, L, B_cell, l, init_B_cell, T] = build_matrix_lifting_shrinking(W, card, init_ML, init_CL, parent_B_cell)

    % W: gram matrix n x n
    % card: cardinality matrix k x k
    % init_ML: matrix of must-link constraints (n_ml x 2)
    % init_CL: matrix of cannot-link constraints (n_cl x 2)
    % parent_Bcell: valid inequalities for the parent node
    % return problem data sdpnal+ format
    
    original_n = size(W, 1);
    k = size(card, 1);
    
    % transformation for shrinking
    [T, ML_graph] = build_T(original_n, init_ML);
    % plot(ML_graph)
    n = size(T, 1); % new size
    CL = update_CL(ML_graph, init_CL);
    e = T * ones(original_n, 1);
    % disp('NEW CL')
    % disp(CL)
    W_full = zeros(n+k, n+k);
    W_full(k+1:n+k, k+1:n+k) = T * W * T';
    
    Y_tl = Y_slice_C(n, k);
    [X_rowsum, X_colsum] = Y_slice_X_shrinking(n, k, e);
    [Z_rowsum, Z_diag] = Y_slice_Z_shrinking(n, k, card, e);
    
    [B_cell_X_cl, l_X_cl, A_cell_Z_cl, b_Z_cl] = add_cannot_link_Y(n, k, CL);

    Acell = [Y_tl, X_rowsum, X_colsum, Z_rowsum, Z_diag, ... 
        A_cell_Z_cl];
    b_card = [];
    for i = 1:k
        for j = i:k
            b_card = [b_card; card(i, j)];
        end
    end
    b = [b_card; ones(n, 1); diag(card); ones(n, 1); zeros(n, 1); ... 
        b_Z_cl];
    
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
    
    
    % inherit cuts and change indices
    init_B_cell = inherit_cuts_matrix_lifting(ML_graph, parent_B_cell, n, k);

end