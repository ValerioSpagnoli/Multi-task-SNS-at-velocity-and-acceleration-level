
function q_ddot_PS = pseudoinverse_solution(m, n, J, J_dot, x_ddot, bounds, q_dot, verbose)

    bounds_min_acceleration = bounds{3}(1,:);
    bounds_max_acceleration = bounds{3}(2,:);


    % number of task
    l = length(m);

    P = {};
    q_ddot = {};
    
    P_0 = eye(n);
    q_ddot_0 = zeros(n,1);
    
    for k=1:l
        if k==1
            P_km1 = P_0;
            q_ddot_km1 = q_ddot_0;
        else
            P_km1 = P{k-1};
            q_ddot_km1 = q_ddot{k-1};            
        end

        J_k = J{k};
        J_dot_k = J_dot{k};
        x_ddot_k = x_ddot{k};
        P{k} = P_km1 - pinv(J_k*P_km1)*J_k*P_km1;

        q_ddot_k = q_ddot_km1 + pinv(J_k*P_km1) * (x_ddot_k - J_dot_k*q_dot - J_k*q_ddot_km1); 

        s = 1;
        s_new = 1;
        for i=1:n
            if q_ddot_k(i) > bounds_max_acceleration(i)
                s_new = bounds_max_acceleration(i)/q_ddot_k(i);
            elseif q_ddot_k(i) < bounds_min_acceleration(i)
                s_new = bounds_min_acceleration(i)/q_ddot_k(i);
            end

            if s_new < s
                s = s_new;
            end
        end

        q_ddot_k = s*q_ddot_k;

        
    end
   
    q_ddot_PS = q_ddot_k;