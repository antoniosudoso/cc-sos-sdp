function [B_triangle, n_ineq] = separate_triangle_Y(Y, k, eps, max_triangle_ineq, triangle_perc)
    
    rng(1727);
    
    sizeY = size(Y, 1); % (n+k)*(n+k)
    n = sizeY - k;
    n_idx = randperm(n);
    
    B_triangle = cell(1, max_triangle_ineq);
    counter_triangle = 1;
    
    violations = zeros(max_triangle_ineq, 1);
    
    stop = false;
    
    % TRIANGLE INEQUALITIES
    
    for i=n_idx
        for j=i+1:n
            for t=j+1:n
                 if (i ~= t)
                     viol1 = 0.5*Y(i+k,j+k) + 0.5*Y(j+k,i+k) + 0.5*Y(i+k,t+k) + 0.5*Y(t+k,i+k) - Y(i+k,i+k) - 0.5*Y(j+k,t+k) - 0.5*Y(t+k,j+k);
                     viol2 = 0.5*Y(i+k,j+k) + 0.5*Y(j+k,i+k) + 0.5*Y(j+k,t+k) + 0.5*Y(t+k,j+k) - Y(j+k,j+k) - 0.5*Y(i+k,t+k) - 0.5*Y(t+k,i+k);
                     viol3 = 0.5*Y(i+k,t+k) + 0.5*Y(t+k,i+k) + 0.5*Y(j+k,t+k) + 0.5*Y(t+k,j+k) - Y(t+k,t+k) - 0.5*Y(i+k,j+k) - 0.5*Y(j+k,i+k);
                     if viol1 >= eps
                        violations(counter_triangle) = viol1;
                        B_triangle{counter_triangle} = sparse([i+k, j+k, i+k, t+k, i+k, j+k, t+k], [j+k, i+k, t+k, i+k, i+k, t+k, j+k], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+k, n+k);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol2 >= eps
                        violations(counter_triangle) = viol2;
                        B_triangle{counter_triangle} = sparse([i+k, j+k, j+k, t+k, j+k, i+k, t+k], [j+k, i+k, t+k, j+k, j+k, t+k, i+k], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+k, n+k);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol3 >= eps
                        violations(counter_triangle) = viol3;
                        B_triangle{counter_triangle} = sparse([i+k, t+k, j+k, t+k, t+k, i+k, j+k], [t+k, i+k, t+k, j+k, t+k, j+k, i+k], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+k, n+k);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     
                     if counter_triangle - 1 == max_triangle_ineq 
                         stop = true;
                         break;
                     end   
                 end
            end
            if stop == true
                break;
            end
        end
        if stop == true
            break;
        end
    end
    
    n_ineq = counter_triangle - 1;
    max_ineq = max_triangle_ineq * triangle_perc;
    if n_ineq <= max_ineq
        B_triangle = B_triangle(1:n_ineq);
    else    
        violations = violations(1:n_ineq);
        [~, id_sorted] = sort(violations, 'descend');
        B_triangle = B_triangle(id_sorted);
        n_ineq = floor(max_ineq);
        B_triangle = B_triangle(1:n_ineq);
        
        clear id_sorted
    end
    
    clear n_idx violations 

end