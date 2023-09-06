% Parameters:
%  > m = {m_1, ..., m_k}: cell array with dimension of each task 
%  > J = {J_1, ..., J_k}: cell array with Jacobian matrix, one for each task 
%  > x_dot = {x_dot_1, ..., x_dot_k}: cell array with desired task velocity
%  > bounds_max = [bounds_max_1, ..., bounds_max_n]: maximum bounds for each joint
%  > bounds_min = [bounds_min_1, ..., bounds_min_n]: minimum bounds for each joint
%  > verbose = <true, false>: enable prints
% Return: q_dot_SNS: solution of multitask SNS algorithm 

function q_ddot_SNS = SNS_acceleration_multitask(n, m, J, J_dot, x_ddot, bounds, q, q_dot, T, verbose)
    
    bounds_min_position = bounds{1}(1,:);
    bounds_max_position = bounds{1}(2,:);

    bounds_min_velocity = bounds{2}(1,:);
    bounds_max_velocity = bounds{2}(2,:);

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);

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

        W_star_k = eye(n);
        P_bar_star_k = eye(n);
        q_ddot_star_N_k = zeros(n,1);
        s_star_k = 0;

        W_k = eye(n);
        q_ddot_N_k = zeros(n,1);
        s_k=1;
        
        if k==1
            P_km1 = P_0;
            q_ddot_bar_km1 = q_ddot_bar_0;
        else
            P_km1 = round( P{k-1} ,4);
            q_ddot_bar_km1 = round( q_ddot_bar{k-1} ,4);
        end 
        P_bar_k = P_km1;
        
        J_k = round( J{k}, 4);
        J_dot_k = round( J_dot{k}, 4);
        x_ddot_k = round( x_ddot{k}, 4);
        m_k = m{k};
        
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
                if while_loop>n-m_k
                    break; 
                end
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
            q_ddot_bar_k = round( q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) ,4);
            
            % shaping velocity bounds
            bounds_Q_ddot_min = zeros(7,1);
            bounds_Q_ddot_max = zeros(7,1);
            for i=1:n

                bounds_Q_ddot_min_1 = round( (2*(bounds_min_position(i)-q(i)-q_dot(i)*T))/T^2 ,4);
                bounds_Q_ddot_min_2 = round( -(bounds_max_velocity(i)+q_dot(i))/T ,4);
                bounds_Q_ddot_min_3 = round( -bounds_max_acceleration(i), 4);
                
                % if max([bounds_Q_ddot_min_1, bounds_Q_ddot_min_2, bounds_Q_ddot_min_3]) < 0 && q_dot(i) < (bounds_min_position(i)-q(i))/T
                %    bounds_Q_ddot_min_2 = -bounds_max_velocity(i)/T - (bounds_min_position(i)-q(i))/T^2;
                % end

                bounds_Q_ddot_min(i) = max([bounds_Q_ddot_min_1, bounds_Q_ddot_min_2, bounds_Q_ddot_min_3]);


                bounds_Q_ddot_max_1 = round( (2*(bounds_max_position(i)-q(i)-q_dot(i)*T))/T^2 ,4);
                bounds_Q_ddot_max_2 = round( (bounds_max_velocity(i)-q_dot(i))/T ,4);
                bounds_Q_ddot_max_3 = round( bounds_max_acceleration(i) ,4);

                % if min([bounds_Q_ddot_max_1, bounds_Q_ddot_max_2, bounds_Q_ddot_max_3]) > 0 && q_dot(i) > (bounds_min_position(i)-q(i))/T
                %    bounds_Q_ddot_max_2 = bounds_max_velocity(i)/T - (bounds_min_position(i)-q(i))/T^2;
                % end

                bounds_Q_ddot_max(i) = min([bounds_Q_ddot_max_1, bounds_Q_ddot_max_2, bounds_Q_ddot_max_3]);                
               
            end

            if verbose
            fprintf('J_k = \n');disp(J_k);
            fprintf('P_bar_k = \n');disp(P_bar_k);
            fprintf('pinv(J_k*P_bar_k) = \n');disp(pinv_Jk_x_Pbark);
            fprintf('q_ddot_bar_km1 = ');disp(q_ddot_bar_km1');
            fprintf('q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
            fprintf('q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
            fprintf('q_ddot_bar_k   = ');disp(q_ddot_bar_k');
            end
        
            % check if q_ddot_bar_k has at least one component out of bounds
            for i=1:n
                if verbose
                fprintf('q_ddot_bar_k(i) = %f, bounds_Q_ddot_max(i) = %f, bounds_Q_ddot_min(i) = %f\n', q_ddot_bar_k(i), bounds_Q_ddot_max(i), bounds_Q_ddot_min(i));
                end
                if q_ddot_bar_k(i) < bounds_Q_ddot_min(i) || q_ddot_bar_k(i) > bounds_Q_ddot_max(i)
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

                a = round( pinv_Jk_x_Pbark * x_ddot_k ,4);
                b = round( q_ddot_bar_k - a ,4);

                if verbose
                fprintf('a = ');disp(a')
                fprintf('b = ');disp(b')
                end

                % compute the task scaling factor and the most critical joint
                [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_ddot_max, bounds_Q_ddot_min, verbose); 
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
                    q_ddot_star_N_k = q_ddot_N_k;
                    P_bar_star_k = P_bar_k;
 
                    if verbose
                    fprintf('..........................................................\n')
                    fprintf('task_scaling_factor > s_star_k:\n');
                    fprintf('    W_star_k = \n');disp(W_star_k);
                    fprintf('    P_bar_star_k = \n');disp(P_bar_star_k);
                    fprintf('    s_stark_k = ');disp(s_star_k);
                    fprintf('    q_dot_star_N_k = ');disp(q_ddot_star_N_k');
                    fprintf('..........................................................\n')
                    end
                       
                end

                % saturate the most critical joint
                W_k(j,j) = 0;
                if q_ddot_bar_k(j) > bounds_Q_ddot_max(j)
                    q_ddot_N_k(j) = round( bounds_Q_ddot_max(j)-q_ddot_bar_km1(j) ,4);
                elseif q_ddot_bar_k(j) < bounds_Q_ddot_min(j)
                    q_ddot_N_k(j) = round( bounds_Q_ddot_min(j)-q_ddot_bar_km1(j) ,4);
                end
                
                P_bar_k = round( (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1 ,4);

                if verbose
                fprintf('..........................................................\n')
                fprintf('recomputing stuff for most critical joint:\n');
                fprintf('    W_k = \n');disp(W_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('    q_ddot_N_k = ');disp(q_ddot_N_k');
                fprintf('    J_k*P_bar_k = \n');disp(J_k*P_bar_k);
                fprintf('    rank(J_k*P_bar_k) = ');disp(rank(round(J_k*P_bar_k,4)));
                fprintf('..........................................................\n')
                end

                % if the rank of J_k*P_bar_k is smaller than the dimension
                % of the current task (m_k), scale the remanining
                % components out of bounds using the best scaling factor
                % found so far
                if rank(round(J_k*P_bar_k,4)) < m_k
                    W_k = W_star_k;
                    P_bar_k = P_bar_star_k;
                    q_ddot_N_k = q_ddot_star_N_k;
                    s_k = s_star_k;

                    q_ddot_bar_N_k = round( pinv((eye(n)-W_k)*P_km1)*q_ddot_N_k ,4);
                    q_ddot_tilde_k = round( q_ddot_bar_km1 + q_ddot_bar_N_k ,4);
                    
                    pinv_Jk_x_Pbark = round( pinv(J_k*P_bar_k) ,4);
                    for j = 1:n
                        if P_bar_k(j,j) == 0
                            pinv_Jk_x_Pbark(j,:) = 0;
                        end
                    end
                    q_ddot_bar_k = round( q_ddot_tilde_k + pinv_Jk_x_Pbark*(s_k*x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_tilde_k) ,4);

                    limit_exceeded = false;   

                    if verbose
                    fprintf('..........................................................\n')
                    fprintf('rank(J_k*P_bar_k) < m_k:\n');
                    fprintf('    s_k = ');disp(s_k)
                    fprintf('    W_k = \n');disp(W_k);
                    fprintf('    P_bar_k = \n');disp(P_bar_k);
                    fprintf('    q_ddot_N_k = ');disp(q_ddot_N_k');
                    fprintf('    q_ddot_bar_N_k = ');disp(q_ddot_bar_N_k');
                    fprintf('    q_ddot_tilde_k = ');disp(q_ddot_tilde_k');
                    fprintf('    q_ddot_bar_k = ');disp(q_ddot_bar_k');
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
        
        q_ddot_bar{length(q_ddot_bar)+1} = q_ddot_bar_k;

        P_k = round( P_km1 - pinv(J_k*P_km1) * (J_k*P_km1) ,4);
        P{length(P)+1} = P_k;
        
    end
    
    q_ddot_SNS = q_ddot_bar_k;


% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_max, bounds_min, verbose)

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