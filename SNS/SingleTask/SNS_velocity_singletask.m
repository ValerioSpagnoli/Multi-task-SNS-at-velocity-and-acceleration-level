% SNS Algorithm 1
function q_dot_SNS = SNS_velocity_singletask(m, n, J, x_dot, bounds_max, bounds_min, verbose)
    W = eye(n);
    q_dot_N = zeros(n,1);
    s = 1;
    s_star = 0;    

    limit_exceeded = true;
    W_star = W;
    q_dot_N_star = q_dot_N;

    while limit_exceeded

        limit_exceeded = false;
        q_dot_SNS = round(q_dot_N + pinv(J*W)*((s*x_dot)-(J*q_dot_N)),4);

        if verbose
        fprintf('##########################################################\n');        
        fprintf('q_dot_N: \n');disp(round(q_dot_N,4));
        fprintf('s: %f\n', round(s,4));
        fprintf('J*W = ');disp(round(J*W,4));
        fprintf('q_dot_SNS: ');disp(round(q_dot_SNS,4));
        end
        
        for i=1:n
            if q_dot_SNS(i) > bounds_max(i) || q_dot_SNS(i) < bounds_min(i)
                limit_exceeded = true;
            end
        end

        if verbose
        fprintf('limit_exceeded: %d\n', limit_exceeded);
        fprintf('----------------------------------------------------------\n');
        end

        if limit_exceeded

            % Algorithm 2   
            [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, J, W, x_dot, q_dot_N, bounds_max, bounds_min);

            if verbose
            fprintf('task_scaling_factor: %f\n', round(task_scaling_factor,4));  
            end

            if task_scaling_factor > s_star
                s_star = task_scaling_factor;                
                W_star = W;
                q_dot_N_star = q_dot_N;
            end                            
            
            j = most_critical_joint;

            if verbose
            fprintf('most_critical_joint: %d\n', j);
            end

            W(j,j) = 0;
            if q_dot_SNS(j) > bounds_max(j)
                q_dot_N(j) = bounds_max(j);
            elseif q_dot_SNS(j) < bounds_min(j)
                q_dot_N(j) = bounds_min(j);
            end

            if verbose
            fprintf('q_dot_N = \n');disp(round(q_dot_N,4));
            fprintf('J*W = ');disp(round(J*W,4));       
            fprintf('rank(J*W) = %d\n', rank(J*W))
            fprintf('----------------------------------------------------------\n');
            end

            if rank(J*W) < m

                s = s_star;
                W = W_star;
                q_dot_N = q_dot_N_star;

                if verbose
                fprintf('Scaling ...\n')  
                fprintf('s = %f\n', round(s,4));
                fprintf('W = \n');disp(round(W,4));
                fprintf('q_dot_N = \n');disp(round(q_dot_N,4));
                end

                q_dot_SNS = round(q_dot_N + pinv(J*W)*((s*x_dot)-(J*q_dot_N)),4);
                limit_exceeded = false;
            end
        end
    end

    q_dot_SNS = round(q_dot_SNS, 4);



% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, J, W, x_dot, q_dot_N, bounds_max, bounds_min)

    a = round(pinv(J*W)*x_dot,4);
    b = round(q_dot_N - pinv(J*W)*J*q_dot_N,4);

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n
        Smin = round((bounds_min(i)-b(i))/a(i),4);     
        Smax = round((bounds_max(i)-b(i))/a(i),4);   

        if isinf(Smax) || isnan(Smax)
            Smax = 1e15;
        end
        if isinf(Smin) || isnan(Smin)
            Smin = -1e15;
        end
        
        if Smin > Smax
            S_min(i) = Smax;
            S_max(i) = Smin;
        else
            S_min(i) = Smin;
            S_max(i) = Smax;
        end    
    end

    [s_max, most_critical_joint] = min(S_max);
    s_min = max(S_min);

    if s_min > s_max || s_max < 0 || s_min > 1
        task_scaling_factor = 0;
    else
        task_scaling_factor = min([s_max,1]);
    end
