function [i_idx, j_idx, max_val] = get_branching_pair(PI, card)
    
    if ~iscell(PI)
        n = size(PI, 1);
        Z = PI;
    else
        n_blocks = size(PI, 2);
        n = size(PI{1}, 1);
        Z = zeros(n, n);
        for j=1:n_blocks
            Z = Z + (1/card(j,j)) * PI{j};
        end   
    end
    
    % [i, j, min_val]
    values = zeros(n*(n-1)/2, 3);
    count = 1;
    for i=1:n
        for j=i+1:n
            norm_val = norm(Z(i, :) - Z(j, :))^2;
            first = Z(i, j);
            min_val = min(first, norm_val);
            values(count, :) = [i, j, min_val];
            count = count + 1;
        end
    end
    
    sorted_values_Z = sortrows(values, 3, 'descend');
    i_idx = sorted_values_Z(1, 1);
    j_idx = sorted_values_Z(1, 2);
    max_val = sorted_values_Z(1, 3);
    
    if max_val < 1e-4
        i_idx = -1;
        j_idx = -1;
    end
    
    if j_idx < i_idx
        temp_idx = j_idx;
        j_idx = i_idx;
        i_idx = temp_idx;
    end

end