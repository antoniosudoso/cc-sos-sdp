function [B_clique, n_ineq, viol_vec, C_list] = separate_clique_Y(Y, k, max_card, eps, max_ineq)
    
    sizeY = size(Y, 1); % (n+k)*(n+k)
    n = sizeY - k;
    X = Y(k+1:n+k, k+1:n+k);

    [row_i, col_j] = find(X < 1/max_card);
    V = unique(sort([row_i, col_j], 2), 'rows', 'stable');
    n_pairs = size(V, 1);
    %disp(n_pairs)

    B_clique = cell(1, n_pairs);
    counter_clique = 1;
    viol_vec = zeros(n_pairs, 1);
    C_list = zeros(n_pairs, k+1); % list of violated cliques
    rhs = 1/max_card;
    
    % Greedy Separation of Clique Inequalities
    for h=randperm(n_pairs)
            C = zeros(k, 1);
            C(1) = V(h, 1);
            C(2) = V(h, 2);
            C_size = 2;
            while C_size <= k
                temp_sum = inf;
                temp_u = inf;
                for v=1:n
                    sum = 0;
                    for i=1:C_size
                        if i ~= v
                            sum = sum + X(C(i), v);
                        else
                            sum = inf;
                            break;
                        end
                    end
                    if sum < temp_sum
                        temp_sum = sum;
                        temp_u = v;
                    end
                end
                C_size = C_size + 1;
                C(C_size) = temp_u;  
            end 
            
            lhs = 0;
            row_idx = zeros(C_size*(C_size-1)/2, 1);
            col_idx = zeros(C_size*(C_size-1)/2, 1);
            count = 1;
            for i=1:C_size
                for z=i+1:C_size
                    lhs = lhs + X(C(i), C(z));
                    row_idx(count) = C(i) + k;
                    col_idx(count) = C(z) + k;
                    count = count + 1;
                end
            end
            
            size_idx = size(row_idx, 1);
            viol = lhs - rhs;

            if viol < -eps
                B_clique{counter_clique} = sparse([row_idx; col_idx], [col_idx; row_idx], 0.5*ones(2 * size_idx, 1), n+k, n+k);
                viol_vec(counter_clique) = viol;
                C_list(counter_clique, :) = C';
                counter_clique = counter_clique + 1;
            end
    end
    
    n_ineq = counter_clique - 1;
    C_list = C_list(1:n_ineq, :);
    C_list = unique(sort(C_list, 2), 'rows', 'stable'); % get unique cliques (<= n_ineq)
    
    if n_ineq <= max_ineq
        B_clique = B_clique(1:n_ineq);
        viol_vec = viol_vec(1:n_ineq);
    else
        viol_vec = viol_vec(1:n_ineq);
        [~, id_sorted] = sort(viol_vec, 'ascend');
        B_clique = B_clique(id_sorted);
        n_ineq = floor(max_ineq);
        B_clique = B_clique(1:n_ineq);
    end
    
end