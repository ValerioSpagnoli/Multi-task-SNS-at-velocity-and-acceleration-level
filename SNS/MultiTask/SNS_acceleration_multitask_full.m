function q_ddot_SNS = SNS_acceleration_multitask_full(n, m, J, J_dot, task_d, bounds, q, q_dot, T, verbose)

    % [bounds_Q_ddot_min, bounds_Q_ddot_max] = shaping_acceleration_bounds(n, bounds, q, q_dot, T);
    bounds_Q_ddot_min = bounds{3}(1,:);
    bounds_Q_ddot_max = bounds{3}(2,:);

    l = length(m);

    P_0 = eye(n);
    q_ddot_bar_0 = zeros(n,1);
    P = {};
    q_ddot_bar = {};
    
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
            q_ddot_bar_km1 = round( q_ddot_bar{k-1} ,4);
        end        

        J_k = round(J{k},4);
        J_dot_k = round(J_dot{k},4);
        task_d_k = round(task_d{k},4);
        m_k = m{k};

        if size(task_d_k,1) == n % CS type task
            q_ddot_bar_k = algorithm_4(n, task_d_k, q_ddot_bar_km1, P_km1, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose);
        else                     % normal task
            [q_ddot_bar_k, ret] = algorithm_3(n, m_k, task_d_k, J_k, J_dot_k, q_ddot_bar_km1, P_km1, q_dot, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose);
        end

        q_ddot_bar{length(q_ddot_bar)+1} = q_ddot_bar_k;

        P_k = round( P_km1 - pinv(J_k*P_km1) * (J_k*P_km1) ,4);
        P{length(P)+1} = P_k;

        if ret
            fprintf('End of SNS multitask algorithm.\n')
            break;
        end
    end
    
    q_ddot_SNS = q_ddot_bar_k;
end


function [q_ddot_bar_k, ret] = algorithm_3(n, m_k, x_ddot_k, J_k, J_dot_k, q_ddot_bar_km1, P_km1, q_dot, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose)

    ret=false;

    W_star_k = eye(n);
    q_ddot_star_N_k = zeros(n,1);
    P_bar_star_k = eye(n);
    s_star_k = 0;

    W_k = eye(n);
    q_ddot_N_k = zeros(n,1);
    s_k=1;
    P_bar_k = P_km1;
    
    if verbose
    fprintf('P_bar_k = P_km1 = \n');disp(P_bar_k);
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

        q_ddot_bar_N_k = round( pinv((eye(n)-W_k)*P_km1)*q_ddot_N_k ,4);
        q_ddot_tilde_k = round( q_ddot_bar_km1 + q_ddot_bar_N_k ,4);
        
        pinv_Jk_x_Pbark = round( pinv(J_k*P_bar_k) ,4);
        for j = 1:n
            if P_bar_k(j,j) == 0
                pinv_Jk_x_Pbark(j,:) = 0;
            end
        end
        q_ddot_bar_k = round(q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) ,4);
        
        if verbose
        fprintf('pinv(J_k*P_bar_k) = \n');disp(pinv_Jk_x_Pbark);
        fprintf('J_dot = \n');disp(J_dot_k);
        fprintf('x_ddot = ');disp(x_ddot_k');
        fprintf('q_dot = '); disp(q_dot');
        fprintf('q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
        fprintf('q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
        fprintf('q_ddot_bar_k = ');disp(q_ddot_bar_k');
        end
    
        for i=1:n
            if verbose
            fprintf('\nq_ddot_bar_k(i) = %f, bounds_Q_ddot_max(i) = %f, bounds_Q_ddot_min(i) = %f\n', q_ddot_bar_k(i), bounds_Q_ddot_max(i), bounds_Q_ddot_min(i));
            end
            if q_ddot_bar_k(i) > bounds_Q_ddot_max(i) || q_ddot_bar_k(i) < bounds_Q_ddot_min(i)
                limit_exceeded = true;
                if verbose
                fprintf(' --> limit_exceeded on joint %d\n',i);
                end
            end
            if verbose
                fprintf('\n');
            end
        end

        if limit_exceeded
            if verbose              
            fprintf('----------------------------------------------------------\n')
            fprintf('limit_exceeded: %d\n\n', limit_exceeded);                
            end

            a = round( pinv_Jk_x_Pbark * x_ddot_k ,4); % accelerazione desiderata per i joint non saturati
            % b = round( q_ddot_bar_k - a ,4);           % differenza tra accelerazione attuale e accelerazione desiderata
            b = round( q_ddot_N_k - pinv_Jk_x_Pbark*(J_dot_k*q_dot+J_k *q_ddot_N_k) ,4);

            if verbose
            fprintf('..........................................................\n')
            end
            [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_ddot_max, bounds_Q_ddot_min, verbose); 
            j = most_critical_joint;

            if verbose                
            fprintf('..........................................................\n')
            fprintf('task_scaling_factor: %f\n', task_scaling_factor);
            fprintf('most_critical_joint: %d\n', j);
            fprintf('..........................................................\n')
            end

            if task_scaling_factor > s_star_k
                s_star_k = task_scaling_factor;
                W_star_k = W_k;
                q_ddot_star_N_k = q_ddot_N_k;
                P_bar_star_k = P_bar_k;

                if verbose
                fprintf('..........................................................\n')
                fprintf('task_scaling_factor > s_star_k:\n');
                fprintf('    W_star_k = \n');disp(W_star_k);
                fprintf('    P_bar_star_k = \n');disp(P_bar_star_k);
                fprintf('    s_star_k = ');disp(s_star_k);
                fprintf('    q_ddot_star_N_k = ');disp(q_ddot_star_N_k');
                fprintf('..........................................................\n')
                end
            end

            W_k(j,j) = 0;
            if q_ddot_bar_k(j) > bounds_Q_ddot_max(j)
                q_ddot_N_k(j) = round( bounds_Q_ddot_max(j)-q_ddot_bar_km1(j) ,4);
            elseif q_ddot_bar_k(j) < bounds_Q_ddot_min(j)
                q_ddot_N_k(j) = round( bounds_Q_ddot_min(j)-q_ddot_bar_km1(j) ,4);
            end
            
            P_bar_k = (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1;

            if verbose
            fprintf('..........................................................\n')
            fprintf('recomputing stuff for most critical joint:\n');
            fprintf('    W_k = \n');disp(W_k);
            fprintf('    P_bar_k = \n');disp(P_bar_k);
            fprintf('    q_ddot_N_k = ');disp(q_ddot_N_k');
            fprintf('..........................................................\n')
            end

            if rank(J_k*P_bar_k) < m_k
                q_ddot_bar_k_pre = round(q_ddot_tilde_k + pinv_Jk_x_Pbark*(x_ddot_k -J_dot_k*q_dot - J_k*q_ddot_tilde_k) ,4);

                s_k = s_star_k;
                W_k = W_star_k;
                q_ddot_N_k = q_ddot_star_N_k;
                P_bar_k = P_bar_star_k;

                q_ddot_bar_N_k = round( pinv((eye(n)-W_k)*P_km1)*q_ddot_N_k ,4);
                q_ddot_tilde_k = round( q_ddot_bar_km1 + q_ddot_bar_N_k ,4);
                
                pinv_Jk_x_Pbark = round(pinv(J_k*P_bar_k),4);
                for j = 1:n
                    if P_bar_k(j,j) == 0
                        pinv_Jk_x_Pbark(j,:) = 0;
                    end
                end
                q_ddot_bar_k = round(q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) ,4);
                
                limit_exceeded = false;       
                ret = true;

                if verbose
                fprintf('..........................................................\n')
                fprintf('rank(J_k*P_bar_k < m_k):\n');
                fprintf('    s_k = ');disp(s_k)
                fprintf('    W_k = \n');disp(W_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('    q_ddot_N_k = ');disp(q_ddot_N_k');
                fprintf('    q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
                fprintf('    q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
                fprintf('    J_dot_k*q_dot = ');disp((J_dot_k*q_dot)');
                fprintf('    J_k*q_ddot_tilde_k = ');disp((J_k*q_ddot_tilde_k)');
                fprintf('    - J_dot_k*q_dot - J_k*q_ddot_tilde_k = ');disp((- J_dot_k*q_dot - J_k*q_ddot_tilde_k)');
                fprintf('    q_ddot_bar_k_pre = ');disp(q_ddot_bar_k_pre');
                fprintf('    q_ddot_bar_k = ');disp(q_ddot_bar_k');
                fprintf('..........................................................\n')
                end
            else
                ret = false;
            end % if rank(J_k*P_bar_k) < m_k

            if verbose
            fprintf('----------------------------------------------------------\n')
            end

        end % if limit_exceeded

        if verbose
        fprintf('**********************************************************\n')
        end
        while_loop = while_loop+1;

    end % while limit_exceeded
end



function q_ddot_SNS = algorithm_4(n, q_ddot_cs, q_ddot_bar_km1, P_km1, bounds_Q_ddot_min, bounds_Q_ddot_max, verbose)
    W_cs = eye(n);
    
    for i=1:n
        if q_ddot_bar_km1(i) == bounds_Q_ddot_min(i) || q_ddot_bar_km1(i) == bounds_Q_ddot_max(i)
            W_cs(i,i) = 0;
        end
    end

    P_bar_cs = (eye(n) - pinv((eye(n)-W_cs)*P_km1))*P_km1;
    a = P_bar_cs*q_ddot_cs;
    b = q_ddot_bar_km1;

    [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_ddot_max, bounds_Q_ddot_min, verbose);
    s_cs = task_scaling_factor;
    q_ddot_SNS = q_ddot_bar_km1 + s_cs*P_bar_cs*q_ddot_cs; 
end




% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_max, bounds_min, verbose)

    % a = pinv_Jk_x_Pbark*(x_ddot_k-J_dot_k*q_dot)  accelerazione desiderata per i joint non saturati
    % b = q_ddot_bar_k - a                          differenza tra accelerazione attuale e accelerazione desiderata

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n
        Smax = round( (bounds_max(i)-b(i))/a(i) ,4);  % bound - b = new_bound 
        Smin = round( (bounds_min(i)-b(i))/a(i) ,4);    

        if isinf(Smax) || isnan(Smax)
            Smax = +Inf;
        end
        if isinf(Smin) || isnan(Smin)
            Smin = -Inf;
        end

        if verbose
        fprintf('bounds_max(%d) = ',i);disp(bounds_max(i));
        fprintf('bounds_min(%d) = ',i);disp(bounds_min(i));
        fprintf('a(%d) = ',i);disp(a(i));
        fprintf('b(%d) = ',i);disp(b(i));
        fprintf('Smax = (bounds_max(%d)-b(%d))/a(%d) = ',i,i,i);disp(Smax);
        fprintf('Smin = (bounds_min(%d)-b(%d))/a(%d) = ',i,i,i);disp(Smin);
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
        task_scaling_factor = s_max; %min([s_max,1]);
    end
end

%% Compute acceleration bounds 
function [bounds_Q_ddot_min, bounds_Q_ddot_max] = shaping_acceleration_bounds(n, bounds, q, q_dot, T)

    bounds_min_position = bounds{1}(1,:);
    bounds_max_position = bounds{1}(2,:);

    bounds_min_velocity = bounds{2}(1,:);
    bounds_max_velocity = bounds{2}(2,:);

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);

    bounds_Q_ddot_min = zeros(1,7);
    bounds_Q_ddot_max = zeros(1,7);

    for i=1:n

        Q_ddot_min_1 = 2*(bounds_min_position(i)-q(i)-q_dot(i)*T)/T^2;
        Q_ddot_min_2 = -(bounds_max_velocity(i)+q_dot(i))/T;
        Q_ddot_min_3 = bounds_min_acceleration(i);
        
        bounds_Q_ddot_min(i) = round( max([Q_ddot_min_1, Q_ddot_min_2, Q_ddot_min_3]) ,4);

        Q_ddot_max_1 = 2*(bounds_max_position(i)-q(i)-q_dot(i)*T)/T^2;
        Q_ddot_max_2 = (bounds_max_velocity(i)-q_dot(i))/T;
        Q_ddot_max_3 = bounds_max_acceleration(i);
        
        bounds_Q_ddot_max(i) = round( min([Q_ddot_max_1, Q_ddot_max_2, Q_ddot_max_3]) ,4);        
    end    
end