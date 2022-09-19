function [B_pair, n_ineq] = separate_pair_Y(Y, k, eps, max_pair_ineq, pair_perc)

    rng(1727);
    
    sizeY = size(Y, 1); % (n+k)*(n+k)
    n = sizeY - k;
    n_idx = randperm(n);
    
    B_pair = cell(1, max_pair_ineq);
    counter_pair = 1;
    
    violations = zeros(max_pair_ineq, 1);
    
    stop = false;
    
    % PAIR INEQUALITIES
    
    for i=n_idx
        for j=i+1:n
            viol1 = 0.5*Y(i+k, j+k) + 0.5*Y(j+k, i+k) - Y(i+k, i+k);
            viol2 = 0.5*Y(i+k, j+k) + 0.5*Y(j+k, i+k) - Y(j+k, j+k);
            if viol1 >= eps
                violations(counter_pair) = viol1;
                B_pair{counter_pair} = sparse([i+k ; j+k ; i+k], [j+k; i+k; i+k], [-0.5; -0.5; 1], n+k, n+k);
                counter_pair = counter_pair + 1;
            end
            if viol2 >= eps
                violations(counter_pair) = viol2;
                B_pair{counter_pair} = sparse([i+k ; j+k ; j+k], [j+k; i+k; j+k], [-0.5; -0.5; 1], n+k, n+k);
                counter_pair = counter_pair + 1;
            end
            if counter_pair - 1 == max_pair_ineq
                stop = true;
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    n_ineq = counter_pair - 1;
    max_ineq = max_pair_ineq * pair_perc;
    if n_ineq <= max_ineq
        B_pair = B_pair(1:n_ineq);
    else
        violations = violations(1:n_ineq);
        [~, id_sorted] = sort(violations, 'descend');
        B_pair = B_pair(id_sorted);
        n_ineq = floor(max_ineq);
        B_pair = B_pair(1:n_ineq);
        
        clear id_sorted
    end
    
    clear violations
    
end
