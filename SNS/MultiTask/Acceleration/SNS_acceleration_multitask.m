function q_ddot_SNS = SNS_acceleration_multitask(n, m, J, J_dot, task_ddot, bounds, q, q_dot, T, verbose)
    
    [bounds_Q_ddot_min, bounds_Q_ddot_max] = shaping_joint_acceleration_bounds(n, bounds, q, q_dot, T);

    l = length(m);

    P_0 = eye(n);
    q_ddot_bar_0 = zeros(n,1);
    P = {};
    q_ddot_bar = {};

    saturated_joints = zeros(1,n);
        
    for k=1:l
        if verbose
        fprintf('##########################################################\n')
        fprintf('task number = %d\n', k);
        fprintf('##########################################################\n')
        end
        
        if k==1
            P_km1 = P_0;
            q_ddot_bar_km1 = q_ddot_bar_0;
        else
            P_km1 = P{k-1};
            q_ddot_bar_km1 = q_ddot_bar{k-1};            
        end 
        P_bar_k = P_km1;
        
        J_k = J{k};
        J_dot_k = J_dot{k};
        task_ddot_k = task_ddot{k};            

        m_k = m{k};

        if m_k < n
            if(norm(task_ddot_k)>0)
                [q_ddot_bar_k, saturated_joints] = algorithm_3(n, m_k, task_ddot_k, J_k, J_dot_k, q_ddot_bar_km1, P_km1, q_dot, bounds_Q_ddot_min, bounds_Q_ddot_max, saturated_joints, verbose);            
            end
        else                        
            q_ddot_bar_k = algorithm_4(n, task_ddot_k, q_ddot_bar_km1, P_km1, bounds_Q_ddot_min, bounds_Q_ddot_max, saturated_joints, verbose);
        end              

        q_ddot_bar{length(q_ddot_bar)+1} = q_ddot_bar_k;
                  
        P_k = round( P_km1 - pinv(J_k*P_km1) * (J_k*P_km1) ,5);
        P{length(P)+1} = P_k;
        
    end
    
    q_ddot_SNS = q_ddot_bar_k;
end



%% SNS Algorithm 3: compute new acceleration from cartesian task
function [q_ddot_bar_k, saturated_joints] = algorithm_3(n, m_k, x_ddot_k, J_k, J_dot_k, q_ddot_bar_km1, P_km1, q_dot, bounds_Q_ddot_min, bounds_Q_ddot_max, saturated_joints, verbose)
    
W_star_k = eye(n);
    P_bar_star_k = eye(n);
    q_ddot_star_N_k = zeros(n,1);
    s_star_k = 0;
    most_critical_joint_star = -1;
    saturated_joints_star = saturated_joints;
    last_saturated_joint = -1;

    W_k = eye(n);
    q_ddot_N_k = zeros(n,1);
    s_k=1;
    P_bar_k = P_km1;    

    if verbose
        fprintf('P_bar_k = P_km1 = \n');disp(P_km1);
        fprintf('q_ddot_bar_km1 = ');disp(q_ddot_bar_km1');
    end

    while_loop = 1;
    limit_exceeded = true;
    while limit_exceeded

        if verbose
        fprintf('**********************************************************\n')
        fprintf('while loop %d\n\n', while_loop);   

        end

        limit_exceeded = false;
               
        q_ddot_bar_N_k = pinv((eye(n)-W_k)*P_km1) * q_ddot_N_k;
        q_ddot_tilde_k = q_ddot_bar_km1 + q_ddot_bar_N_k;

        pinv_Jk_x_Pbark = pinv(J_k*P_bar_k);
        for j = 1:n
            if saturated_joints(j) == 1
                pinv_Jk_x_Pbark(j,:) = 0;
            end
        end
        q_ddot_bar_k = q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k);
        
        if verbose
        fprintf('saturated_joints = ');disp(saturated_joints);
        fprintf('s_k = ');disp(s_k);
        fprintf('x_ddot_k = ');disp(x_ddot_k');
        fprintf('q_dot = ');disp(q_dot');
        fprintf('\n');
        fprintf('J_k = \n');disp(J_k);
        fprintf('J_dot_k = \n');disp(J_dot_k);
        fprintf('P_km1 = \n');disp(P_km1);
        fprintf('P_bar_k = \n');disp(P_bar_k);
        fprintf('pinv(J_k*P_bar_k) = \n');disp(pinv_Jk_x_Pbark);
        fprintf('\n');
        fprintf('q_ddot_N_k  = ');disp(q_ddot_N_k');
        fprintf('q_ddot_bar_km1 = ');disp(q_ddot_bar_km1');
        fprintf('q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
        fprintf('q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
        fprintf('\n');
        fprintf('s_k*x_ddot_k = ');disp((s_k*x_ddot_k)')
        fprintf('J_dot_k*q_dot = ');disp((J_dot_k*q_dot)');
        fprintf('J_k*q_ddot_tilde_k = ');disp((J_k*q_ddot_tilde_k)');
        fprintf('pinv(J_k*P_bar_k)*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) = ');disp((pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k))');
        fprintf('\n');
        fprintf('q_ddot_bar_k   = ');disp(q_ddot_bar_k');           
        end
    
        % check if q_ddot_bar_k has at least one component out of bounds
        for i=1:n
            if verbose
            fprintf('q_ddot_bar_k(i) = %f, bounds_Q_ddot_max(i) = %f, bounds_Q_ddot_min(i) = %f\n', q_ddot_bar_k(i), bounds_Q_ddot_max(i), bounds_Q_ddot_min(i));
            end
            if q_ddot_bar_k(i) < round(bounds_Q_ddot_min(i),4) || q_ddot_bar_k(i) > round(bounds_Q_ddot_max(i),4)
                if verbose
                fprintf(' --> limit_exceeded on joint %d\n', i);
                end
                limit_exceeded = true;
            end
        end

        % if q_ddot_bar_k has at least one component out of bounds
        if limit_exceeded

            if verbose              
            fprintf('----------------------------------------------------------\n')
            fprintf('limit_exceeded: %d\n\n', limit_exceeded);                
            end

            a = pinv_Jk_x_Pbark * x_ddot_k;            
            b = q_ddot_tilde_k - pinv_Jk_x_Pbark * J_dot_k * q_dot - pinv_Jk_x_Pbark * J_k * q_ddot_tilde_k;
                                    
            if verbose
            fprintf('a = ');disp(a')
            fprintf('b = ');disp(b')
            end

            % compute the task scaling factor and the most critical joint
            [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose);            

            if verbose                
            fprintf('..........................................................\n')
            fprintf('task_scaling_factor: %f\n', task_scaling_factor);
            fprintf('most_critical_joint: %d\n', most_critical_joint);
            fprintf('..........................................................\n')
            end

            % if the computed task scaling factor is greater than the
            % current best task scaling factor, save the current
            % configuration
            if task_scaling_factor > s_star_k
                saved_best_config = true;
                s_star_k = task_scaling_factor;
                most_critical_joint_star = most_critical_joint;
                W_star_k = W_k;
                q_ddot_star_N_k = q_ddot_N_k;
                P_bar_star_k = P_bar_k;
                saturated_joints_star = saturated_joints;

                if verbose
                fprintf('..........................................................\n')
                fprintf('task_scaling_factor > s_star_k:\n');
                fprintf('    W_star_k = \n');disp(W_star_k);
                fprintf('    P_bar_star_k = \n');disp(P_bar_star_k);
                fprintf('    s_stark_k = ');disp(s_star_k);
                fprintf('    q_ddot_star_N_k = ');disp(q_ddot_star_N_k');
                fprintf('    saturated_joints_star = ');disp(saturated_joints_star);
                fprintf('..........................................................\n')
                end    
            end

            j = most_critical_joint;
            last_saturated_joint = j;
            saturated_joints(j) = 1;

            % saturate the most critical joint
            W_k(j,j) = 0;
            if q_ddot_bar_k(j) > bounds_Q_ddot_max(j)
                q_ddot_N_k(j) = bounds_Q_ddot_max(j)-q_ddot_bar_km1(j);
            elseif q_ddot_bar_k(j) < bounds_Q_ddot_min(j)
                q_ddot_N_k(j) = bounds_Q_ddot_min(j)-q_ddot_bar_km1(j);
            end            
            
            P_bar_k = round( (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1 ,8);                      

            if verbose
            fprintf('..........................................................\n')
            fprintf('recomputing stuff for most critical joint:\n');
            fprintf('    saturated_joints = ');disp(saturated_joints);
            fprintf('    W_k = \n');disp(W_k);
            fprintf('    P_bar_k = \n');disp(P_bar_k);
            fprintf('    q_ddot_N_k = ');disp(q_ddot_N_k');
            fprintf('    J_k*P_bar_k = \n');disp(J_k*P_bar_k);
            fprintf('    rank(J_k*P_bar_k) = ');disp(rank(J_k*P_bar_k));
            fprintf('..........................................................\n')
            end

            
            % if the rank of J_k*P_bar_k is smaller than the dimension
            % of the current task (m_k), scale the remanining
            % components out of bounds using the best scaling factor
            % found so far
            if rank(J_k*P_bar_k) < m_k
                s_k = s_star_k;
                most_critical_joint = most_critical_joint_star;
                W_k = W_star_k;
                P_bar_k = P_bar_star_k;
                q_ddot_N_k = q_ddot_star_N_k;
                saturated_joints = saturated_joints_star;  

                q_ddot_bar_N_k = pinv((eye(n)-W_k)*P_km1) * q_ddot_N_k;
                q_ddot_tilde_k = q_ddot_bar_km1 + q_ddot_bar_N_k;                
                pinv_Jk_x_Pbark = pinv(J_k*P_bar_k);
                for j = 1:n
                    if saturated_joints(j) == 1
                        pinv_Jk_x_Pbark(j,:) = 0;
                    end
                end
                q_ddot_bar_k = q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k);

                if most_critical_joint ~= -1
                    saturated_joints(most_critical_joint) = 1;
                end  

                limit_exceeded = false;   

                if verbose
                fprintf('..........................................................\n')
                fprintf('rank(J_k*P_bar_k) < m_k:\n');
                fprintf('    s_k = ');disp(s_k)
                fprintf('    x_ddot_k = ');disp(x_ddot_k');
                fprintf('    q_dot = ');disp(q_dot');
                fprintf('    W_k = \n');disp(W_k);     
                fprintf('    saturated_joints = ');disp(saturated_joints);
                fprintf('\n');
                fprintf('    J_k = \n');disp(J_k);
                fprintf('    J_dot_k = \n');disp(J_dot_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('    pinv(J_k*P_bar_k) = \n');disp(pinv_Jk_x_Pbark);
                fprintf('\n');
                fprintf('    q_ddot_N_k  = ');disp(q_ddot_N_k');
                fprintf('    q_ddot_bar_km1 = ');disp(q_ddot_bar_km1');
                fprintf('    q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
                fprintf('    q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
                fprintf('\n');
                fprintf('    s_k*x_ddot_k = ');disp((s_k*x_ddot_k)')
                fprintf('    J_dot_k*q_dot = ');disp((J_dot_k*q_dot)');
                fprintf('    J_k*q_ddot_tilde_k = ');disp((J_k*q_ddot_tilde_k)');
                fprintf('    pinv(J_k*P_bar_k)*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) = ');disp((pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k))');
                fprintf('\n');
                fprintf('    q_ddot_bar_k   = ');disp(q_ddot_bar_k'); 
                fprintf('..........................................................\n')
                end                
            end

            if verbose
            fprintf('----------------------------------------------------------\n')
            end

        end
        if verbose
        fprintf('**********************************************************\n')
        end
        while_loop = while_loop+1;
    end
end



%% SNS Algorithm 4: compute new acceleration from a configuration space task
function q_ddot_bar_k = algorithm_4(n, q_ddot_cs, q_ddot_bar_km1, P_km1, bounds_Q_ddot_min, bounds_Q_ddot_max, saturated_joints, verbose)

    W_cs = eye(n);

    for j=1:n
        if q_ddot_bar_km1(j) == bounds_Q_ddot_min(j) || q_ddot_bar_km1(j) == bounds_Q_ddot_max(j)
            W_cs(j,j) = 0;
        end
    end
    
    P_bar_cs = (eye(n) - pinv((eye(n)-W_cs)*P_km1))*P_km1;
   
    a = P_bar_cs*q_ddot_cs;
    b = q_ddot_bar_km1;

    [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose);    
    s_cs = task_scaling_factor;
    q_ddot_bar_k = q_ddot_bar_km1 + s_cs*P_bar_cs*q_ddot_cs;   
end



%% SNS algorithm 2: get scaling factor
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_min, bounds_max, verbose)

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n
        Smin = (bounds_min(i)-b(i))/a(i);
        Smax = (bounds_max(i)-b(i))/a(i);        

        if isinf(Smax) || isnan(Smax)
            Smax = +Inf;
        end
        if isinf(Smin) || isnan(Smin)
            Smin = -Inf;
        end

        if verbose
        fprintf('bounds_max(i) = ');disp(bounds_max(i));
        fprintf('bounds_min(i) = ');disp(bounds_min(i));
        fprintf('b(i) = ');disp(b(i));
        fprintf('a(i) = ');disp(a(i));

        fprintf('Smax = (bounds_max(i)-b(i))/a(i) = ');disp(Smax);
        fprintf('Smin = (bounds_min(i)-b(i))/a(i) = ');disp(Smin);
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

    if verbose
    fprintf('S_max = ');disp(S_max);
    fprintf('S_min = ');disp(S_min);

    fprintf('s_max = ');disp(s_max);
    fprintf('s_min = ');disp(s_min);
    fprintf('most_critical_joint = ');disp(most_critical_joint);
    end

    if s_min > s_max || s_max < 0 || s_min > 1
        task_scaling_factor = 0;
    else
        task_scaling_factor = min([s_max,1]);
    end
end



%% Shaping joint acceleration bounds
function [bounds_Q_ddot_min, bounds_Q_ddot_max] = shaping_joint_acceleration_bounds(n, bounds, q, q_dot, T)

    bounds_min_position = bounds{1}(1,:);
    bounds_max_position = bounds{1}(2,:);

    bounds_min_velocity = bounds{2}(1,:);
    bounds_max_velocity = bounds{2}(2,:);

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);

    bounds_Q_ddot_min = zeros(7,1);
    bounds_Q_ddot_max = zeros(7,1);
    for i=1:n        
        bounds_Q_ddot_min_1 = min([(2*(bounds_min_position(i)-q(i)-q_dot(i)*T))/T^2, 0]);
        bounds_Q_ddot_min_2 = min([-(bounds_max_velocity(i)+q_dot(i))/T,0]);
        bounds_Q_ddot_min_3 = -bounds_max_acceleration(i);

        bounds_Q_ddot_min(i) = max([bounds_Q_ddot_min_1, bounds_Q_ddot_min_2, bounds_Q_ddot_min_3]);

        bounds_Q_ddot_max_1 = max([(2*(bounds_max_position(i)-q(i)-q_dot(i)*T))/T^2,0]);
        bounds_Q_ddot_max_2 = max([(bounds_max_velocity(i)-q_dot(i))/T,0]);
        bounds_Q_ddot_max_3 = bounds_max_acceleration(i);

        bounds_Q_ddot_max(i) = min([bounds_Q_ddot_max_1, bounds_Q_ddot_max_2, bounds_Q_ddot_max_3]);          
    end
end