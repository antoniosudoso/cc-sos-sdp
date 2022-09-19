function [blk, At, C, b, L] = build_matrix_lifting(data, card)

    % data: data matrix n x d
    % card: cardinality matrix k x k
    % return problem data sdpnal+ format

    n = size(data, 1);
    k = size(card, 1);
    
    D = 0.5 * squareform(pdist(data, 'squaredeuclidean'));

    W_full = zeros(n+k, n+k);
    W_full(k+1:n+k, k+1:n+k) = -D;
    
    Y_tl = Y_slice_C(n, k);
    [X_rowsum, X_colsum] = Y_slice_X(n, k);
    [Z_rowsum, Z_diag] = Y_slice_Z(n, k, card);

    Acell = [Y_tl, X_rowsum, X_colsum, Z_rowsum, Z_diag];
    b_card = [];
    for i = 1:k
        for j = i:k
            b_card = [b_card; card(i, j)];
        end
    end
    b = [b_card; ones(n, 1); diag(card); ones(n, 1); zeros(n, 1)];
    blk = cell(1);
    blk{1,1} = 's';
    blk{1,2} = n+k;

    At = svec(blk, Acell, 1);

    C = cell(1);
    C{1} = -sparse(W_full);

    L = cell(1);
    L{1} = 0; % we want Y >= 0

end