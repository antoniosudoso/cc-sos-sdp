function [B_pair, n_ineq] = separate_pair_P(Y, eps, max_pair_ineq, pair_perc)
    
    rng(1727);
    
    sizeY = size(Y, 1); % (n+1)*(n+1)
    n = sizeY - 1;
    n_idx = randperm(n);
    
    B_pair = cell(1, max_pair_ineq);
    counter_pair = 1;
    
    violations = zeros(max_pair_ineq, 1);
    
    stop = false;
    
    % PAIR INEQUALITIES
    
    for i=n_idx
        for j=i+1:n
            viol1 = 0.5*Y(i+1, j+1) + 0.5*Y(j+1, i+1) - Y(i+1, i+1);
            viol2 = 0.5*Y(i+1, j+1) + 0.5*Y(j+1, i+1) - Y(j+1, j+1);
            %viol3 = Y(i+1, i+1) + Y(j+1, j+1) - 0.5*Y(i+1, j+1) - 0.5*Y(j+1, i+1);
            if viol1 >= eps
                violations(counter_pair) = viol1;
                B_pair{counter_pair} = sparse([i+1; j+1; i+1], [j+1; i+1; i+1], [-0.5; -0.5; 1], n+1, n+1);
                counter_pair = counter_pair + 1;
            end
            if viol2 >= eps
                violations(counter_pair) = viol2;
                B_pair{counter_pair} = sparse([i+1; j+1; j+1], [j+1; i+1; j+1], [-0.5; -0.5; 1], n+1, n+1);
                counter_pair = counter_pair + 1;
            end
            %if viol3 >= eps
            %    violations(counter_pair) = viol3;
            %    B_pair{counter_pair} = sparse([i+1; j+1; i+1; j+1], [i+1; j+1; j+1; i+1], [-1; -1; 0.5; 0.5], n+1, n+1);
            %    counter_pair = counter_pair + 1;
            %end
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