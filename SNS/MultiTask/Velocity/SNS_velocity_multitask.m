% Parameters:
%  > m = {m_1, ..., m_k}: cell array with dimension of each task 
%  > J = {J_1, ..., J_k}: cell array with Jacobian matrix, one for each task 
%  > x_dot = {x_dot_1, ..., x_dot_k}: cell array with desired task velocity
%  > bounds_max = [bounds_max_1, ..., bounds_max_n]: maximum bounds for each joint
%  > bounds_min = [bounds_min_1, ..., bounds_min_n]: minimum bounds for each joint
%  > verbose = <true, false>: enable prints
% Return: q_dot_SNS: solution of multitask SNS algorithm 


function q_dot_SNS = SNS_velocity_multitask(n, m, J, task_dot, bounds, q, T, verbose)
    

    if verbose
    fprintf('##########################################################\n')
    end
    
    [bounds_Q_dot_min, bounds_Q_dot_max] = shaping_joint_velocity_bounds(n, bounds, q, T, verbose);
    
    if verbose
    fprintf('##########################################################\n')
    end

    l = length(m);

    P_0 = eye(n);
    q_dot_bar_0 = zeros(n,1);
    P = {};
    q_dot_bar = {};
        
    for k=1:l

        if verbose
        fprintf('##########################################################\n')
        fprintf('task number = %d\n', k);
        fprintf('##########################################################\n')
        end
        
        if k==1
            P_km1 = P_0;
            q_dot_bar_km1 = q_dot_bar_0;
        else
            P_km1 = round( P{k-1} ,4);
            q_dot_bar_km1 = round( q_dot_bar{k-1} ,4);
        end 
        P_bar_k = P_km1;
        
        J_k = round( J{k}, 4);
        task_dot_k = round( task_dot{k}, 4);
        m_k = m{k};

        if m_k < n
            if norm(task_dot_k)>0
                q_dot_bar_k = algorithm_3(n, m_k, task_dot_k, J_k, q_dot_bar_km1, P_km1, bounds_Q_dot_min, bounds_Q_dot_max, verbose);
            end
        else            
            q_dot_bar_k = algorithm_4(n, task_dot_k, q_dot_bar_km1, P_km1, bounds_Q_dot_min, bounds_Q_dot_max, verbose);
        end
        
        q_dot_bar{length(q_dot_bar)+1} = q_dot_bar_k;
        
        P_k = round( P_km1 - pinv(J_k*P_km1) * (J_k*P_km1) ,4);
        P{length(P)+1} = P_k;
        
    end
    
    q_dot_SNS = q_dot_bar_k;
end



% SNS Algorithm 3
function q_dot_bar_k = algorithm_3(n, m_k, x_dot_k, J_k, q_dot_bar_km1, P_km1, bounds_Q_dot_min, bounds_Q_dot_max, verbose)
    W_star_k = eye(n);
    P_bar_star_k = eye(n);
    q_dot_star_N_k = zeros(n,1);
    s_star_k = 0;

    W_k = eye(n);
    q_dot_N_k = zeros(n,1);
    s_k=1;
    P_bar_k = P_km1;

    if verbose
        fprintf('P_bar_k = P_km1 = \n');disp(P_km1);
        fprintf('q_dot_bar_km1 = ');disp(q_dot_bar_km1');
    end

    while_loop = 1;
    limit_exceeded = true;
    while limit_exceeded

        if verbose
        fprintf('**********************************************************\n')
        fprintf('while loop %d\n\n', while_loop);   

        end

        limit_exceeded = false;
               
        pinv_IW_x_Pkm1 = round( pinv((eye(n)-W_k)*P_km1) , 4);
        for j = 1:n
            if q_dot_N_k(j) ~= 0
                pinv_IW_x_Pkm1(j,j) = 1;
            end
        end
        q_dot_bar_N_k = round( pinv_IW_x_Pkm1 * q_dot_N_k ,4);

        q_dot_tilde_k = round( q_dot_bar_km1 + q_dot_bar_N_k ,4);

        pinv_Jk_x_Pbark = round( pinv(J_k*P_bar_k) ,4);
        for j = 1:n
            if P_bar_k(j,j) == 0
                pinv_Jk_x_Pbark(j,:) = 0;
            end
        end
        q_dot_bar_k = round( q_dot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_dot_k - J_k*q_dot_tilde_k) ,4);
        

        if verbose
        fprintf('s_k = ');disp(s_k);
        fprintf('x_dot_k = ');disp(x_dot_k');
        fprintf('J_k = \n');disp(J_k);
        fprintf('P_bar_k = \n');disp(P_bar_k);
        fprintf('pinv(J_k*P_bar_k) = \n');disp(pinv_Jk_x_Pbark);
        fprintf('pinv((eye(n)-W_k)*P_km1) = \n');disp(pinv((eye(n)-W_k)*P_km1));
        fprintf('q_dot_N_k  = ');disp(q_dot_N_k');
        fprintf('q_dot_bar_km1 = ');disp(q_dot_bar_km1');        
        fprintf('q_dot_bar_N_k = ');disp(q_dot_bar_N_k');
        fprintf('q_dot_tilde_k = ');disp(q_dot_tilde_k');
        fprintf('q_dot_bar_k   = ');disp(q_dot_bar_k');        
        end
    
        % check if q_dot_bar_k has at least one component out of bounds
        for i=1:n
            if verbose
            fprintf('q_dot_bar_k(i) = %f, bounds_Q_dot_max(i) = %f, bounds_Q_dot_min(i) = %f\n', q_dot_bar_k(i), bounds_Q_dot_max(i), bounds_Q_dot_min(i));
            end
            if q_dot_bar_k(i) < bounds_Q_dot_min(i) || q_dot_bar_k(i) > bounds_Q_dot_max(i)
                if verbose
                fprintf(' --> limit_exceeded on joint %d\n', i);
                end
                limit_exceeded = true;
            end
        end

        % if q_dot_bar_k has at least one component out of bounds
        if limit_exceeded

            if verbose              
            fprintf('----------------------------------------------------------\n')
            fprintf('limit_exceeded: %d\n\n', limit_exceeded);                
            end

            a = round( pinv_Jk_x_Pbark * x_dot_k ,4);
            % b = round( q_dot_bar_k - a ,4);
            b = q_dot_N_k - pinv_Jk_x_Pbark * J_k * q_dot_N_k;

            if verbose
            fprintf('a = ');disp(a')
            fprintf('b = ');disp(b')
            end

            % compute the task scaling factor and the most critical joint
            [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_dot_min, bounds_Q_dot_max, verbose); 
            j = most_critical_joint;

            if verbose                
            fprintf('..........................................................\n')
            fprintf('task_scaling_factor: %f\n', task_scaling_factor);
            fprintf('most_critical_joint: %d\n', j);
            fprintf('..........................................................\n')
            end


            % if the computed task scaling factor is greater than the
            % current best task scaling factor, save the current
            % configuration
            if task_scaling_factor > s_star_k
                s_star_k = task_scaling_factor;
                W_star_k = W_k;
                q_dot_star_N_k = q_dot_N_k;
                P_bar_star_k = P_bar_k;

                if verbose
                fprintf('..........................................................\n')
                fprintf('task_scaling_factor > s_star_k:\n');
                fprintf('    W_star_k = \n');disp(W_star_k);
                fprintf('    P_bar_star_k = \n');disp(P_bar_star_k);
                fprintf('    s_star_k = ');disp(s_star_k);
                fprintf('    q_dot_star_N_k = ');disp(q_dot_star_N_k');
                fprintf('..........................................................\n')
                end    
            end


            % saturate the most critical joint
            W_k(j,j) = 0;
            if q_dot_bar_k(j) > bounds_Q_dot_max(j)
                q_dot_N_k(j) = bounds_Q_dot_max(j)-q_dot_bar_km1(j);
            elseif q_dot_bar_k(j) < bounds_Q_dot_min(j)
                q_dot_N_k(j) = bounds_Q_dot_min(j)-q_dot_bar_km1(j);
            end
            
            P_bar_k = round( (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1 ,4);

            if verbose
            fprintf('..........................................................\n')
            fprintf('recomputing stuff for most critical joint:\n');
            fprintf('    W_k = \n');disp(W_k);
            fprintf('    P_bar_k = \n');disp(P_bar_k);
            fprintf('    q_dot_N_k = ');disp(q_dot_N_k');
            fprintf('    J_k*P_bar_k = \n');disp(J_k*P_bar_k);
            fprintf('    rank(J_k*P_bar_k) = ');disp(rank(J_k*P_bar_k));
            fprintf('..........................................................\n')
            end

            
            % if the rank of J_k*P_bar_k is smaller than the dimension
            % of the current task (m_k), scale the remanining
            % components out of bounds using the best scaling factor
            % found so far
            if rank(J_k*P_bar_k) < m_k
                W_k = W_star_k;
                P_bar_k = P_bar_star_k;
                q_dot_N_k = q_dot_star_N_k;
                s_k = s_star_k;


                pinv_IW_x_Pkm1 = round( pinv((eye(n)-W_k)*P_km1) , 4);
                for j = 1:n
                    if q_dot_N_k(j) ~= 0
                        pinv_IW_x_Pkm1(j,j) = 1;
                    end
                end
                q_dot_bar_N_k = round( pinv_IW_x_Pkm1 * q_dot_N_k ,4);

                q_dot_tilde_k = round( q_dot_bar_km1 + q_dot_bar_N_k ,4);
                
                pinv_Jk_x_Pbark = round( pinv(J_k*P_bar_k) ,4);
                for j = 1:n
                    if P_bar_k(j,j) == 0
                        pinv_Jk_x_Pbark(j,:) = 0;
                    end
                end
                q_dot_bar_k = round( q_dot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_dot_k - J_k*q_dot_tilde_k) ,4);

                limit_exceeded = false;   

                if verbose
                fprintf('..........................................................\n')
                fprintf('rank(J_k*P_bar_k) < m_k:\n');
                fprintf('    s_k = ');disp(s_k)
                fprintf('    W_k = \n');disp(W_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('    q_dot_N_k = ');disp(q_dot_N_k');
                fprintf('    q_dot_bar_N_k = ');disp(q_dot_bar_N_k');
                fprintf('    q_dot_tilde_k = ');disp(q_dot_tilde_k');
                fprintf('    q_dot_bar_k = ');disp(q_dot_bar_k');
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



% SNS Algorithm 4
function q_dot_bar_k = algorithm_4(n, q_dot_cs, q_dot_bar_km1, P_km1, bounds_Q_dot_min, bounds_Q_dot_max, verbose)

    if verbose
        fprintf('P_bar_k = P_km1 = \n');disp(P_km1);
        fprintf('q_ddot_bar_km1 = ');disp(q_dot_bar_km1');
    end

    W_cs = eye(n);

    for i=1:n
        if q_dot_bar_km1(i) == bounds_Q_dot_min(i) || q_dot_bar_km1(i) == bounds_Q_dot_max(i)
            W_cs(i,i) = 0;
        end
    end

    P_bar_cs = (eye(n) - pinv((eye(n)-W_cs)*P_km1))*P_km1;
    a = P_bar_cs*q_dot_cs;
    b = q_dot_bar_km1;

    [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_dot_min, bounds_Q_dot_max, verbose);
    s_cs = task_scaling_factor;
    q_dot_bar_k = q_dot_bar_km1 + s_cs*P_bar_cs*q_dot_cs;   
end



% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_min, bounds_max, verbose)

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n

        Smin = round( (bounds_min(i)-b(i))/a(i) ,4);
        Smax = round( (bounds_max(i)-b(i))/a(i) ,4);   

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



% SNS shaping joint velocities bounds
function [bounds_Q_dot_min, bounds_Q_dot_max] = shaping_joint_velocity_bounds(n, bounds, q, T, verbose)
    bounds_min_position = bounds{1}(1,:);
    bounds_max_position = bounds{1}(2,:);

    bounds_min_velocity = bounds{2}(1,:);
    bounds_max_velocity = bounds{2}(2,:);

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);

    % shaping velocity bounds
    bounds_Q_dot_min = zeros(7,1);
    bounds_Q_dot_max = zeros(7,1);

    if verbose
    fprintf('Shaping joint velocities bounds: ');
    end

    for i=1:n

        if verbose
        fprintf(' - Bound %d\n\n',i);

        fprintf('   - bounds_min_position(i) = ');disp(bounds_min_position(i))
        fprintf('   - bounds_max_position(i) = ');disp(bounds_max_position(i))
        fprintf('   - bounds_min_velocity(i) = ');disp(bounds_min_velocity(i))
        fprintf('   - bounds_max_velocity(i) = ');disp(bounds_max_velocity(i))
        fprintf('   - bounds_min_acceleration(i) = ');disp(bounds_min_acceleration(i))
        fprintf('   - bounds_max_acceleration(i) = ');disp(bounds_max_acceleration(i))
        fprintf('   - q(i) = ');disp(q(i));

        fprintf('   - Computations ... \n')

        fprintf('   - bounds_Q_dot_min: \n')
        fprintf('     1. (bounds_min_position(i)-q(i))/T = ');disp((bounds_min_position(i)-q(i))/T);
        fprintf('     2. bounds_min_velocity(i) = ');disp(bounds_min_velocity(i));
        fprintf('     3. -sqrt(2*bounds_max_acceleration(i)*(q(i)-bounds_min_position(i))) = ');disp(-sqrt(2*bounds_max_acceleration(i)*(q(i)-bounds_min_position(i))));
        fprintf('     - bounds_Q_dot_min(i) = max([1, 2, 3]) = ');disp(max([(bounds_min_position(i)-q(i))/T, bounds_min_velocity(i), -sqrt(2*bounds_max_acceleration(i)*(q(i)-bounds_min_position(i)))]));
        
        fprintf('   - bounds_Q_dot_max: \n')
        fprintf('     1. (bounds_max_position(i)-q(i))/T = ');disp((bounds_max_position(i)-q(i))/T);
        fprintf('     2. bounds_max_velocity(i) = ');disp(bounds_max_velocity(i));
        fprintf('     3. sqrt(2*bounds_max_acceleration(i)*(bounds_max_position(i)-q(i))) = ');disp(sqrt(2*bounds_max_acceleration(i)*(bounds_max_position(i)-q(i))));
        fprintf('     - bounds_Q_dot_max(i) = min([1, 2, 3]) = ');disp(min([(bounds_max_position(i)-q(i))/T, bounds_max_velocity(i), sqrt(2*bounds_max_acceleration(i)*(bounds_max_position(i)-q(i)))]));

        fprintf('----------------------------------------------------------\n')
        end

        bounds_Q_dot_min(i) = max([(bounds_min_position(i)-q(i))/T, bounds_min_velocity(i), -sqrt(2*bounds_max_acceleration(i)*(q(i)-bounds_min_position(i)))]);
        bounds_Q_dot_max(i) = min([(bounds_max_position(i)-q(i))/T, bounds_max_velocity(i), sqrt(2*bounds_max_acceleration(i)*(bounds_max_position(i)-q(i)))]);
    end
end