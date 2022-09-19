function [ml_violated, cl_violated] = check_pairwise_lb(X, PI, ML, CL, tol)


    n_ml = size(ML, 1);
    n_cl = size(CL, 1);
    k = size(X, 2);
    n_blocks = k;
    
    if ~iscell(PI)
        temp = PI;
        PI = cell(1);
        PI{1} = temp;
        n_blocks = 1;
    end
    
    ml_violated = 0;
    
    for c=1:n_ml
        i = ML(c, 1);
        j = ML(c, 2);
        fprintf('ML %d %d\n', i, j);
        if norm(X(i, :) - X(j, :)) > tol
            ml_violated = 1;
            break
        end
        for t=1:n_blocks
            viol = norm(PI{t}(i, :) - PI{t}(j, :));
            if viol > tol
                ml_violated = 1;
                break
            end
        end 
        if ml_violated == 1
            break
        end
    end
    
    cl_violated = 0;

    for c=1:n_cl
        i = CL(c, 1);
        j = CL(c, 2);
        fprintf('CL %d %d\n', i, j);
        if sum(abs(X(i, :) + X(j, :)) - 1 <= tol) < k
            cl_violated = 1;
            break
        end
        for t=1:n_blocks
            if abs(PI{t}(i, j)) > tol
                cl_violated = 1;
                break
            end
        end
        if cl_violated == 1
            break
        end
    end

end