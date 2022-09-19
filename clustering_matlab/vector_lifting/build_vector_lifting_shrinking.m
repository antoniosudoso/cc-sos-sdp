function [blk, At, C, b, L, B_cell, l, n_ineq_block, init_B_cell, T] = build_vector_lifting_shrinking(W, card, init_ML, init_CL, parent_B_cell)

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
    D = T * W * T';
    
    n_cl = size(CL, 1);
    
    blk = [];
    L = [];
    C = [];

    At = cell(1);
    A = cell(1);
    b = [];
    
    B_cell = cell(1, k);
    %Bt = cell(1);
    l = [];
    n_ineq_block = zeros(k, 1);
    
    for j = 1:k

        blk{j, 1} = 's';
        blk{j, 2} = n+1;
        L{j, 1} = 0;
        D_temp = zeros(n+1, n+1);
        D_temp(2:n+1, 2:n+1) = D;
        C{j, 1} = -(1 / card(j,j)) * D_temp;
        
        [B_cell_pi_cl, l_pi_cl, A_cell_PI_cl, b_PI_cl] = add_cannot_link_P(n, CL);

        [A_one, A_sum, A_cross_sum, A_rowsum, A_diag] = Y_slice_shrinking(n, card(j,j), e);
        A{j} = [repmat(Y_zeros_cell(2*(n+1)+size(b_PI_cl, 1), n+1), 1, (j-1)), ...
            A_one, A_sum, A_rowsum, A_diag, A_cell_PI_cl, ...
            repmat(Y_zeros_cell(2*(n+1)+size(b_PI_cl, 1), n+1), 1, (k-j)), ...
            A_cross_sum];
        
        At(j) = svec(blk(j, :), A{j}, 1);

        b_temp = [1; card(j,j); zeros(n, 1); zeros(n, 1); b_PI_cl];
        b = [b; b_temp];
        
        B_cell{j} = [B_cell{j}, B_cell_pi_cl];
        n_ineq_block(j) = n_ineq_block(j) + n_cl;
        l = [l; l_pi_cl];

    end
    
    b = [b; ones(n, 1)];
    
    % inherit cuts and change indices
    init_B_cell = inherit_cuts_vector_lifting(ML_graph, parent_B_cell, n, k);

end