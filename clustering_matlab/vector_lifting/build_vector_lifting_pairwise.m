function [blk, At, C, b, L, B_cell, l, n_ineq_block] = build_vector_lifting_pairwise(data, card, ML, CL)

    % data: data matrix n x d
    % card: cardinality matrix k x k
    % return problem data sdpnal+ format
    
    n = size(data, 1);
    k = size(card, 1);
    
    n_ml = size(ML, 1);
    n_cl = size(CL, 1);

    D = 0.5 * squareform(pdist(data, 'squaredeuclidean'));
    
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
        C{j, 1} = (1 / card(j,j)) * D_temp;
        
        [A_cell_pi_ml, b_pi_ml, A_cell_PI_ml, b_PI_ml] = add_must_link_P(n, ML);
        [B_cell_pi_cl, l_pi_cl, A_cell_PI_cl, b_PI_cl] = add_cannot_link_P(n, CL);

        [A_one, A_sum, A_cross_sum, A_rowsum, A_diag] = Y_slice(n, card(j,j));
        A{j} = [repmat(Y_zeros_cell(2*(n+1)+size(b_pi_ml, 1)+size(b_PI_ml, 1)+size(b_PI_cl, 1), n+1), 1, (j-1)), ...
            A_one, A_sum, A_rowsum, A_diag, A_cell_pi_ml, A_cell_PI_ml, A_cell_PI_cl, ...
            repmat(Y_zeros_cell(2*(n+1)+size(b_pi_ml, 1)+size(b_PI_ml, 1)+size(b_PI_cl, 1), n+1), 1, (k-j)), ...
            A_cross_sum];
        
        At(j) = svec(blk(j, :), A{j}, 1);

        b_temp = [1; card(j,j); zeros(n, 1); zeros(n, 1); b_pi_ml; b_PI_ml; b_PI_cl];
        b = [b; b_temp];
        
        B_cell{j} = [B_cell{j}, B_cell_pi_cl];
        n_ineq_block(j) = n_ineq_block(j) + n_cl;
        l = [l; l_pi_cl];

    end
    
    %for j=1:k
    %    B = [Y_zeros_cell(sum(n_ineq_block(1:(j-1))), n+1), ...
    %        B_cell{j}, ...
    %        Y_zeros_cell(sum(n_ineq_block(j+1:k)), n+1)];
    %    Bt(j) = svec(blk(j, :), B, 1);
    %end
    
    b = [b; ones(n, 1)];

end