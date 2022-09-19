function [A_one, A_sum, A_cross_sum, A_rowsum, A_diag] = Y_slice(n, card_j)

    % top-left entry is 1
    A_one = cell(1, 1);
    A_one{1} = sparse(1, 1, 1, n+1, n+1);
    %disp(full(A_one{1}))
    
    % sum vector x
    A_sum = cell(1, 1);
    A_sum{1} = sparse([2:n+1, repelem(1, n)], [repelem(1, n), 2:n+1], 0.5, n+1, n+1);
    %disp(full(A_sum{1}))
    
    % sum across x_j
    A_cross_sum = cell(1, n);
    
    % Pi row sum and diag constraints
    A_rowsum = cell(1, n);
    A_diag = cell(1, n);
    
    c = 1;
    for i = 2:n+1
        coeff = [repelem(0.5, n), repelem(0.5, n), - 0.5 * card_j, - 0.5 * card_j];
        A_rowsum{c} = sparse([repelem(i, n), 2:n+1, i, 1], [2:n+1, repelem(i, n), 1, i], coeff, n+1, n+1);
        %disp(full(A_rowsum{c}))
        coeff = [-0.5, -0.5, 1];
        A_diag{c} = sparse([i, 1, i], [1, i, i], coeff, n+1, n+1);
        %disp(full(A_diag{c}))
        A_cross_sum{c} = sparse([i, 1], [1, i], 0.5, n+1, n+1);
        %disp(full(A_cross_sum{c}))
        c = c + 1;
    end
    
end