% Parameters:
%  > m = {m_1, ..., m_k}: cell array with dimension of each task 
%  > J = {J_1, ..., J_k}: cell array with Jacobian matrix, one for each task 
%  > x_dot = {x_dot_1, ..., x_dot_k}: cell array with desired task velocity
%  > bounds_max = [bounds_max_1, ..., bounds_max_n]: maximum bounds for each joint
%  > bounds_min = [bounds_min_1, ..., bounds_min_n]: minimum bounds for each joint
%  > verbose = <true, false>: enable prints
% Return: q_dot_SNS: solution of multitask SNS algorithm 

function q_dot_SNS = SNS_velocity_multitask(m, n, J, x_dot, bounds, q, T, verbose)
    
    bounds_min_position = bounds{1}(1,:);
    bounds_max_position = bounds{1}(2,:);

    bounds_min_velocity = bounds{2}(1,:);
    bounds_max_velocity = bounds{2}(2,:);

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);

    % number of task
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

        W_star_k = eye(n);
        q_dot_star_N_k = zeros(n,1);
        P_bar_star_k = eye(n);

        W_k = eye(n);
        q_dot_N_k = zeros(n,1);
        s_k=1;
        s_star_k = 0;
        if k==1
            P_km1 = P_0;
            q_dot_bar_km1 = q_dot_bar_0;
        else
            P_km1 = P{k-1};
            q_dot_bar_km1 = round( q_dot_bar{k-1} ,4);
        end 
        P_bar_k = P_km1;
        
        J_k = round(J{k},4);
        x_dot_k = round(x_dot{k},4);
        m_k = m{k};
        q_k = q;
        
        if verbose
        fprintf('W_k = \n');disp(W_k);
        fprintf('q_dot_N_k = \n');disp(q_dot_N_k);
        fprintf('s_k = ');disp(s_k);
        fprintf('s_star_k = ');disp(s_star_k);
        fprintf('P_bar_k = \n');disp(P_bar_k);
        fprintf('P_km1 = \n');disp(P_km1);
        fprintf('q_dot_bar_km1 = \n');disp(q_dot_bar_km1);
        end

        while_loop = 0;
        limit_exceeded = true;
        while limit_exceeded

            if verbose
            fprintf('**********************************************************\n')
            fprintf('while loop %d\n\n', while_loop);            
            end

            limit_exceeded = false;

            q_dot_N_k = round(q_dot_N_k,4);
            
            q_dot_bar_N_k = round( pinv((eye(n)-W_k)*P_km1)*q_dot_N_k ,4);
            q_dot_tilde_k = round( q_dot_bar_km1 + q_dot_bar_N_k ,4);
            q_dot_bar_k = round( q_dot_tilde_k + round(pinv(J_k*P_bar_k),4)*(x_dot_k - J_k*q_dot_tilde_k) ,4);

%             JxP_bar=J_k*P_bar_k;
%             JxP_bar_pinv = round( transpose(JxP_bar)*((JxP_bar*transpose(JxP_bar))^-1) ,4);
%             q_dot_bar_k = round( q_dot_tilde_k + JxP_bar_pinv*(x_dot_k - J_k*q_dot_tilde_k) ,4);
            
            q_dot_bar{length(q_dot_bar)+1} = q_dot_bar_k;

            bounds_Q_dot_min = zeros(7,1);
            bounds_Q_dot_max = zeros(7,1);
            for i=1:n
                bounds_Q_dot_min(i) = max([(bounds_min_position(i)-q_k(i))/T, bounds_min_velocity(i), -sqrt(2*bounds_max_acceleration(i)*(q_k(i)-bounds_min_position(i)))]);
                bounds_Q_dot_max(i) = min([(bounds_max_position(i)-q_k(i))/T, bounds_max_velocity(i), sqrt(2*bounds_max_acceleration(i)*(bounds_max_position(i)-q_k(i)))]);
            end

            if verbose
            fprintf('q_dot_bar_N_k = \n');disp(q_dot_bar_N_k);
            fprintf('q_dot_tilde_k = \n');disp(q_dot_tilde_k);
            fprintf('pinv(J{k}) = \n');disp(pinv(J_k));
            fprintf('pinv(P_bar_k) = \n');disp(pinv(P_bar_k));
            fprintf('J{k}*P_bar_k = \n');disp(J_k*P_bar_k);
            fprintf('pinv(J{k}*P_bar_k) = \n');disp(pinv(round(J_k,4)*round(P_bar_k,4)));
            fprintf('q_dot_bar_k = \n');disp(q_dot_bar_k);
            end
        
            for i=1:n
                if verbose
                fprintf('q_dot_bar_k(i) = %f, bounds_Q_dot_max(i) = %f, bounds_Q_dot_min(i) = %f\n', q_dot_bar_k(i), bounds_Q_dot_max(i), bounds_Q_dot_min(i));
                end
                if q_dot_bar_k(i) > bounds_Q_dot_max(i) || q_dot_bar_k(i) < bounds_Q_dot_min(i)
                    if verbose
                    fprintf(' >> limit_exceeded on joint %d\n', i);
                    end
                    limit_exceeded = true;
                end
            end


            if limit_exceeded

                if verbose              
                fprintf('----------------------------------------------------------\n')
                fprintf('limit_exceeded: %d\n\n', limit_exceeded);                
                end

                a = round( pinv(J_k*P_bar_k)*x_dot_k ,4);
                b = round( q_dot_bar_k - a ,4);

                if verbose
                fprintf('a = \n');disp(a)
                fprintf('b = \n');disp(b)
                end

                [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_Q_dot_max, bounds_Q_dot_min, verbose); 
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
                    q_dot_star_N_k = q_dot_N_k;
                    P_bar_star_k = P_bar_k;
 
                    if verbose
                    fprintf('..........................................................\n')
                    fprintf('task_scaling_factor > s_star_k:\n');
                    fprintf('    s_stark_k = %f\n', s_star_k);
                    fprintf('    W_star_k = \n');disp(W_star_k);
                    fprintf('    q_dot_star_N_k = \n');disp(q_dot_star_N_k);
                    fprintf('    P_bar_star_k = \n');disp(P_bar_star_k);
                    fprintf('..........................................................\n')
                    end
                       
                end


                W_k(j,j) = 0;
                if q_dot_bar_k(j) > bounds_Q_dot_max(j)
                    q_dot_N_k(j) = round( bounds_Q_dot_max(j)-q_dot_bar_km1(j) ,4);
                elseif q_dot_bar_k(j) < bounds_Q_dot_min(j)
                    q_dot_N_k(j) = round( bounds_Q_dot_min(j)-q_dot_bar_km1(j) ,4);
                end
                
                P_bar_k = round( (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1 ,4);

                if verbose
                fprintf('..........................................................\n')
                fprintf('recomputing stuff for most critical joint:\n');
                fprintf('    W_k = \n');disp(W_k);
                fprintf('    q_dot_N_k = \n');disp(q_dot_N_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('..........................................................\n')
                end


                if rank(J_k*P_bar_k) < m_k
                    s_k = s_star_k;
                    W_k = W_star_k;
                    q_dot_N_k = q_dot_star_N_k;
                    P_bar_k = P_bar_star_k;

                    q_dot_bar_N_k = round( pinv((eye(n)-W_k)*P_km1)*q_dot_N_k ,4);
                    q_dot_tilde_k = round( q_dot_bar_km1 + q_dot_bar_N_k ,4);
                    q_dot_bar_k = round( q_dot_tilde_k + round(pinv(J_k*P_bar_k),4)*(s_k*x_dot_k - J_k*q_dot_tilde_k) ,4);

%                     JxP_bar=J_k*P_bar_k;
%                     JxP_bar_pinv = round( transpose(JxP_bar)*((JxP_bar*transpose(JxP_bar))^-1) ,4);
%                     q_dot_bar_k = round( q_dot_tilde_k + JxP_bar_pinv*(s_k*x_dot_k - J_k*q_dot_tilde_k) ,4);

                    limit_exceeded = false;   

                    if verbose
                    fprintf('..........................................................\n')
                    fprintf('rank(J_k*P_bar_k) < m_k:\n');
                    fprintf('    s_k = %f\n', s_k);
                    fprintf('    W_k = \n');disp(W_k);
                    fprintf('    q_dot_N_k = \n');disp(q_dot_N_k);
                    fprintf('    P_bar_k = \n');disp(P_bar_k);
                    fprintf('    q_dot_bar_N_k = \n');disp(q_dot_bar_N_k);
                    fprintf('    q_dot_tilde_k = \n');disp(q_dot_tilde_k);
                    fprintf('    q_dot_bar_k = \n');disp(q_dot_bar_k);
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
        
        P_k = round( P_km1 - round(pinv(J_k*P_km1),4) * (J_k*P_km1) ,4);
        P{length(P)+1} = P_k;
    end
    
    q_dot_SNS = q_dot_bar_k;


% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_max, bounds_min, verbose)

    S_min = zeros(1,n);
    S_max = zeros(1,n);
    
    for i=1:n

        Smax = round( (bounds_max(i)-b(i))/a(i) ,4); 
        Smin = round( (bounds_min(i)-b(i))/a(i) ,4);  

        if isinf(Smax) || isnan(Smax)
            Smax = 1e15;
        end
        if isinf(Smin) || isnan(Smin)
            Smin = -1e15;
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