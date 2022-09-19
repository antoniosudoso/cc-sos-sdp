function [A_one, A_sum, A_cross_sum, A_rowsum, A_diag] = Y_slice_shrinking(m, card_j, e)

    % e = T^s 1_n (new size is m)

    % top-left entry is 1
    A_one = cell(1, 1);
    A_one{1} = sparse(1, 1, 1, m+1, m+1);
    %disp(full(A_one{1}))
    
    % sum vector x
    A_sum = cell(1, 1);
    A_sum{1} = sparse([2:m+1, repelem(1, m)], [repelem(1, m), 2:m+1], [0.5 .* ones(m, 1) .* e; 0.5 .* ones(m, 1) .* e], m+1, m+1);
    %disp(full(A_sum{1}))
    
    % sum across x_j
    A_cross_sum = cell(1, m);
    
    % Pi row sum and diag constraints
    A_rowsum = cell(1, m);
    A_diag = cell(1, m);
    
    c = 1;
    for i = 2:m+1
        coeff = [0.5 .* ones(m, 1) .* e; 0.5 .* ones(m, 1) .* e; - 0.5 * card_j; - 0.5 * card_j];
        A_rowsum{c} = sparse([repelem(i, m), 2:m+1, i, 1], [2:m+1, repelem(i, m), 1, i], coeff, m+1, m+1);
        %disp(full(A_rowsum{c}))
        coeff = [-0.5, -0.5, 1];
        A_diag{c} = sparse([i, 1, i], [1, i, i], coeff, m+1, m+1);
        %disp(full(A_diag{c}))
        A_cross_sum{c} = sparse([i, 1], [1, i], 0.5, m+1, m+1);
        %disp(full(A_cross_sum{c}))
        c = c + 1;
    end
    
end