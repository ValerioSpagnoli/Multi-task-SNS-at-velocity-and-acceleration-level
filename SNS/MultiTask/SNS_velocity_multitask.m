% Parameters:
%  > m = {m_1, ..., m_k}: cell array with dimension of each task 
%  > J = {J_1, ..., J_k}: cell array with Jacobian matrix, one for each task 
%  > x_dot = {x_dot_1, ..., x_dot_k}: cell array with desired task velocity
%  > bounds_max = [bounds_max_1, ..., bounds_max_n]: maximum bounds for each joint
%  > bounds_min = [bounds_min_1, ..., bounds_min_n]: minimum bounds for each joint
%  > verbose = <true, false>: enable prints
% Return: q_dot_SNS: solution of multitask SNS algorithm 

function q_dot_SNS = SNS_velocity_multitask(m, n, J, x_dot, bounds_max, bounds_min, verbose)

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

        W_k = eye(n);
        q_dot_N_k = zeros(n,1);
        s_k=1;
        s_star_k = 0;
        if k==1
            P_km1 = P_0;
            P_bar_k = P_km1;            
            q_dot_bar_km1 = q_dot_bar_0;
        else
            P_bar_k = P{k-1};
            P_km1 = P{k-1};
            q_dot_bar_km1 = q_dot_bar{k-1};
        end        
        
        if verbose
        fprintf('W_k = \n');disp(W_k);
        fprintf('q_dot_N_k = \n');disp(q_dot_N_k);
        fprintf('s_k = ');disp(s_k);
        fprintf('s_star_k = ');disp(s_star_k);
        fprintf('P_bar_k = \n');disp(P_bar_k);
        fprintf('P_km1 = \n');disp(P_km1);
        fprintf('q_dot_bar_km1 = \n');disp(q_dot_bar_km1);
        while_loop = 0;
        end
   
        limit_exceeded = true;
        while limit_exceeded

            if verbose
            fprintf('**********************************************************\n')
            fprintf('while loop %d\n\n', while_loop);            
            end

            limit_exceeded = false;

            q_dot_bar_N_k = pinv((eye(n)-W_k)*P_km1)*q_dot_N_k;
            q_dot_tilde_k = q_dot_bar_km1 + q_dot_bar_N_k;
            q_dot_bar_k = q_dot_tilde_k + pinv(J{k}*P_bar_k)*(s_k*x_dot{k}-J{k}*q_dot_tilde_k);

            q_dot_bar{length(q_dot_bar)+1} = q_dot_bar_k;

            if verbose
            fprintf('q_dot_bar_N_k = \n');disp(q_dot_bar_N_k);
            fprintf('q_dot_tilde_k = \n');disp(q_dot_tilde_k);
            fprintf('q_dot_bar_k = \n');disp(q_dot_bar_k);
            end
        
            for i=1:n
                if q_dot_bar_k(i) > bounds_max(i) || q_dot_bar_k(i) < bounds_min(i)
                    limit_exceeded = true;
                end
            end

            if limit_exceeded

                if verbose              
                fprintf('----------------------------------------------------------\n')
                fprintf('limit_exceeded: %d\n\n', limit_exceeded);                
                end

                a = pinv(J{k}*P_bar_k)*x_dot{k};
                b = q_dot_bar_k - a;
                [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_max, bounds_min); 
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
                if q_dot_bar_k(j) > bounds_max(j)
                    q_dot_N_k(j) = bounds_max(j)-q_dot_bar_km1(j);
                elseif q_dot_bar_k(j) < bounds_min(j)
                    q_dot_N_k(j) = bounds_min(j)-q_dot_bar_km1(j);
                end
                
                P_bar_k = (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1;

                if verbose
                fprintf('..........................................................\n')
                fprintf('recomputing stuff for most critical joint:\n');
                fprintf('    W_k = \n');disp(W_k);
                fprintf('    q_dot_N_k = \n');disp(q_dot_N_k);
                fprintf('    P_bar_k = \n');disp(P_bar_k);
                fprintf('..........................................................\n')
                end

                if rank(J{k}*P_bar_k) < m{k}
                    s_k = s_star_k;
                    W_k = W_star_k;
                    q_dot_N_k = q_dot_star_N_k;
                    P_bar_k = P_bar_star_k;

                    q_dot_bar_N_k = pinv((eye(n)-W_k)*P_km1)*q_dot_N_k;
                    %P_bar_k = (eye(n) - pinv((eye(n)-W_k)*P_km1))*P_km1;
                    q_dot_tilde_k = q_dot_bar_km1 + q_dot_bar_N_k;
                    q_dot_bar_k = q_dot_tilde_k + pinv(J{k}*P_bar_k)*(s_k*x_dot{k}-J{k}*q_dot_tilde_k);

                    limit_exceeded = false;   

                    if verbose
                    fprintf('..........................................................\n')
                    fprintf('rank(J_k*P_bar_k < m_k):\n');
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
            while_loop = while_loop+1;
            fprintf('**********************************************************\n')
            end
            
        end

        P_k = P_km1 - pinv(J{k}*P_km1) * (J{k}*P_km1);
        P{length(P)+1} = P_k;
    end
    
    q_dot_SNS = q_dot_bar_k;


% SNS Algorithm 2
function [task_scaling_factor, most_critical_joint] = getTaskScalingFactor(n, a, b, bounds_max, bounds_min)

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