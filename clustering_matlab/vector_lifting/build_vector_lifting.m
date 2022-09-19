function [blk, At, C, b, L] = build_vector_lifting(data, card)

    % data: data matrix n x d
    % card: cardinality matrix k x k
    % return problem data sdpnal+ format
    
    n = size(data, 1);
    k = size(card, 1);

    D = 0.5 * squareform(pdist(data, 'squaredeuclidean'));

    blk = [];
    L = [];
    C = [];

    At = cell(1);
    A = cell(1);
    b = [];
    for j = 1:k

        blk{j, 1} = 's';
        blk{j, 2} = n+1;
        L{j, 1} = 0;
        D_temp = zeros(n+1, n+1);
        D_temp(2:n+1, 2:n+1) = D;
        C{j, 1} = (1 / card(j,j)) * D_temp;

        [A_one, A_sum, A_cross_sum, A_rowsum, A_diag] = Y_slice(n, card(j,j));
        A{j} = [repmat(Y_zeros_cell(2*(n+1), n+1), 1, (j-1)), ...
            A_one, A_sum, A_rowsum, A_diag, repmat(Y_zeros_cell(2*(n+1), n+1), 1, (k-j)), ...
            A_cross_sum];

        At(j) = svec(blk(j, :), A{j}, 1);

        b_temp = [1; card(j,j); zeros(n, 1); zeros(n, 1)];
        b = [b; b_temp];

    end
    b = [b; ones(n, 1)];

end