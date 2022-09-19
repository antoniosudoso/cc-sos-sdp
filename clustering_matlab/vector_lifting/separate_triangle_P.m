function [B_triangle, n_ineq] = separate_triangle_P(Y, eps, max_triangle_ineq, triangle_perc)
    
    rng(1727);
    
    sizeY = size(Y, 1); % (n+1)*(n+1)
    n = sizeY - 1;
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
                     viol1 = 0.5*Y(i+1,j+1) + 0.5*Y(j+1,i+1) + 0.5*Y(i+1,t+1) + 0.5*Y(t+1,i+1) - Y(i+1,i+1) - 0.5*Y(j+1,t+1) - 0.5*Y(t+1,j+1);
                     viol2 = 0.5*Y(i+1,j+1) + 0.5*Y(j+1,i+1) + 0.5*Y(j+1,t+1) + 0.5*Y(t+1,j+1) - Y(j+1,j+1) - 0.5*Y(i+1,t+1) - 0.5*Y(t+1,i+1);
                     viol3 = 0.5*Y(i+1,t+1) + 0.5*Y(t+1,i+1) + 0.5*Y(j+1,t+1) + 0.5*Y(t+1,j+1) - Y(t+1,t+1) - 0.5*Y(i+1,j+1) - 0.5*Y(j+1,i+1);
                     if viol1 >= eps
                        violations(counter_triangle) = viol1;
                        B_triangle{counter_triangle} = sparse([i+1, j+1, i+1, t+1, i+1, j+1, t+1], ...
                            [j+1, i+1, t+1, i+1, i+1, t+1, j+1], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+1, n+1);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol2 >= eps
                        violations(counter_triangle) = viol2;
                        B_triangle{counter_triangle} = sparse([i+1, j+1, j+1, t+1, j+1, i+1, t+1], ...
                            [j+1, i+1, t+1, j+1, j+1, t+1, i+1], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+1, n+1);
                        counter_triangle = counter_triangle + 1;
                        
                     end
                     if viol3 >= eps
                        violations(counter_triangle) = viol3;
                        B_triangle{counter_triangle} = sparse([i+1, t+1, j+1, t+1, t+1, i+1, j+1], ...
                            [t+1, i+1, t+1, j+1, t+1, j+1, i+1], [-0.5, -0.5, -0.5, -0.5, 1, 0.5, 0.5], n+1, n+1);
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