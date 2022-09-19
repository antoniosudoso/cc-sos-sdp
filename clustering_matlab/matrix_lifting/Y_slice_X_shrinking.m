function [A_rowsum, A_colsum] = Y_slice_X_shrinking(m, k, e)

    % e = T^s 1_n (new size is m)

    A_rowsum = cell(1, m);
    A_colsum = cell(1, k);
    
    c = 1;
    for i = (k+1):(m+k)
        A_rowsum{c} = sparse([repelem(i, k), 1:k], [1:k, repelem(i, k)], 0.5, m+k, m+k);
        %disp(full(A_rowsum{c}))
        c = c + 1;
    end
    
    c = 1;
    for j = 1:k
        A_colsum{c} = sparse([(k+1:m+k)'; repmat(j, m, 1)], [repmat(j, m, 1); (k+1:m+k)'], [0.5 .* ones(m, 1) .* e; 0.5 .* ones(m, 1) .* e], m+k, m+k);
        %disp(full(A_colsum{c}))
        c = c + 1;
    end
    
end