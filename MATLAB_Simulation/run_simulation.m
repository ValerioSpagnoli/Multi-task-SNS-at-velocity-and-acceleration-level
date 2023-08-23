function qs = run_simulation(robot, points, sampling_time)

    ndof = robot.ndof;
    bounds = {robot.bounds_position, robot.bounds_velocity, robot.bounds_acceleration};

    q_0 = robot.q_0;
    q_dot_0 = robot.q_dot_0;
    q_ddot_0 = robot.q_ddot_0;
    ee_position_0 = robot.ee_position_0;
    J_0 = robot.get_J(q_0);
    J_dot_0 = robot.get_J_dot(q_0, q_dot_0);

%     robot_arm = robot{1};
%     bounds = robot{2};
%     ndof = robot{3};
%     q_0 = robot{4};
%     q_dot_0 = zeros(7,1);
%     q_ddot_0 = zeros(7,1);
%     geom_J = geometricJacobian(robot_arm, q_0', 'iiwa_link_ee');
%     J_0 = geom_J(4:6,:);
%     ee_position_0 = tform2trvec(getTransform(robot_arm,q_0','iiwa_link_ee'))';
    
    
    %% initialization

    T = sampling_time;
    kp = 10;
    kd = 0.1;
    
    % current configuration (time h)
    q = q_0;
    q_dot = q_dot_0;
    q_ddot = q_ddot_0;
    ee_position = ee_position_0;
    J = J_0;
    
    % configuration at time h-1
    q_prev = q_0;
    q_dot_prev = q_dot_0;
    q_ddot_prev = q_ddot_0;
    ee_position_prev = ee_position_0;
    J_prev = J_0;
    
    % all configurations during task
    qs = [q_0];

    %% while loop    
    k = 1;
    while true
        
        x = round(points(1:3,k),4);
        m = {length(x)};

        fprintf('k = %d, norm(ee_position - x) = %f\n', k, norm(ee_position-x));
    
        if norm(ee_position - x) < 0.01
            k = k+1;
            x = points(k);
        end
        if k == length(points)+1
            break;
        end

        % compute new configuration
        V = round(kp*norm(x-ee_position) - kd*norm(J_prev*q_dot_prev),4);
        x_dot = round(V * ((x-ee_position)/norm(x-ee_position)),4);
        q_dot_new = SNS_velocity_multitask(m, ndof, {J}, {x_dot}, bounds, q, T, false);
        
        % step forward
        dq = q_dot_new*T;
        q_new = q + dq;

        % set previous configuration to current configuration ...
        q_prev = q;
        q_dot_prev = q_dot;
        J_prev = J;

        % update current cunfiguration to new configuration
        q = q_new;
        q_dot = q_dot_new;
        
        J = robot.get_J(q);
        %geom_J = geometricJacobian(robot_arm, q', 'iiwa_link_ee');
        %J = geom_J(4:6,:);
        
        ee_position = robot.get_ee_position(q);
        %ee_position = tform2trvec(getTransform(robot_arm,q','iiwa_link_ee'))';

        % save new configuration
        qs = [qs, q];
    end
end