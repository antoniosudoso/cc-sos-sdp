function [A_rowsum, A_colsum] = Y_slice_X(n, k)

    A_rowsum = cell(1, n);
    A_colsum = cell(1, k);
    
    c = 1;
    for i = (k+1):(n+k)
        A_rowsum{c} = sparse([repelem(i, k), 1:k], [1:k, repelem(i, k)], 0.5, n+k, n+k);
        %disp(full(A_rowsum{c}))
        c = c + 1;
    end
    
    c = 1;
    for j = 1:k
        A_colsum{c} = sparse([k+1:n+k, repelem(j, n)], [repelem(j, n), k+1:n+k], 0.5, n+k, n+k);
        %disp(full(A_colsum{c}))
        c = c + 1;
    end
    
end