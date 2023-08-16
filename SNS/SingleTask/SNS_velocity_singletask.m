% SNS Algorithm 1
function q_dot_SNS = SNS_velocity_singletask(m, n, J, x_dot, bounds_max, bounds_min)
    W = eye(n);
    q_dot_N = zeros(n,1);
    s = 1;
    s_star = 0;    

    limit_exceeded = true;
    W_star = W;
    q_dot_N_star = q_dot_N;

    while limit_exceeded

        fprintf('##########################################################\n');

        limit_exceeded = false;
        
        fprintf('q_dot_N: \n');disp(round(q_dot_N,4));
        fprintf('s: %f\n', round(s,4));
        fprintf('J*W = ');disp(round(J*W,4));
        
        q_dot_SNS = q_dot_N + pinv(J*W)*((s*x_dot)-(J*q_dot_N));

        fprintf('q_dot_SNS: ');disp(round(q_dot_SNS,4));

        for i=1:n
            if q_dot_SNS(i) > bounds_max(i) || q_dot_SNS(i) < bounds_min(i)
                limit_exceeded = true;
            end
        end

        fprintf('limit_exceeded: %d\n', limit_exceeded);

        fprintf('----------------------------------------------------------\n');

        if limit_exceeded

            % Algorithm 2   
            [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, J, W, x_dot, q_dot_N, bounds_max, bounds_min);
            fprintf('task_scaling_factor: %f\n', round(task_scaling_factor,4));  

            if task_scaling_factor > s_star
                s_star = task_scaling_factor;                
                W_star = W;
                q_dot_N_star = q_dot_N;
            end                            
            
            j = most_critical_joint;
            fprintf('most_critical_joint: %d\n', j);

            W(j,j) = 0;
            if q_dot_SNS(j) > bounds_max(j)
                q_dot_N(j) = bounds_max(j);
            elseif q_dot_SNS(j) < bounds_min(j)
                q_dot_N(j) = bounds_min(j);
            end

            fprintf('q_dot_N = \n');disp(round(q_dot_N,4));
            fprintf('J*W = ');disp(round(J*W,4));       

            fprintf('rank(J*W) = %d\n', rank(J*W))

            fprintf('----------------------------------------------------------\n');

            if rank(J*W) < m
                fprintf('Scaling ...\n')             

                s = s_star;
                W = W_star;
                q_dot_N = q_dot_N_star;

                fprintf('s = %f\n', round(s,4));
                fprintf('W = \n');disp(round(W,4));
                fprintf('q_dot_N = \n');disp(round(q_dot_N,4));

                q_dot_SNS = q_dot_N + pinv(J*W)*((s*x_dot)-(J*q_dot_N));
                limit_exceeded = false;
            end
        end
    end

    q_dot_SNS = round(q_dot_SNS, 4);



% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, J, W, x_dot, q_dot_N, bounds_max, bounds_min)

    a = pinv(J*W)*x_dot;
    b = q_dot_N - pinv(J*W)*J*q_dot_N;

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n
        Smin = (bounds_min(i)-b(i))/a(i);     
        Smax = (bounds_max(i)-b(i))/a(i);         
        
        if Smin > Smax
            temp = Smin;
            S_min = Smax;
            Smax = temp;
        end

        S_min(i) = Smin;
        S_max(i) = Smax;
    end

    [s_max, most_critical_joint] = min(S_max);
    s_min = max(S_min);

    if s_min > s_max || s_max < 0 || s_min > 1
        task_scaling_factor = 0;
    else
        task_scaling_factor = min([s_max,1]);
    end
