function A = Y_slice_C(n, k)

    A = cell(1, k*(k+1)/2);
    c = 1;
    for i = 1:k
        A{c} = sparse(i, i, 1, n+k, n+k);
        %disp(full(A{c}))
        c = c + 1;
        for j = i+1:k
            A{c} = sparse([i, j], [j, i], 0.5, n+k, n+k);
            %disp(full(A{c}))
            c = c + 1;
        end
    end
end