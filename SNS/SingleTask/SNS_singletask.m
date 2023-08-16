function SNS_solution = SNS_singletask(level, m, n, J, J_dot, q_dot, desired_task, bounds_max, bounds_min)
    if strcmp(level, 'velocity')
        SNS_solution = SNS_velocity(m, n, J, desired_task, bounds_max, bounds_min);
    elseif strcmp(level, 'acceleration')
        SNS_solution = SNS_acceleration(m, n, J, J_dot, q_dot, desired_task, bounds_max, bounds_min);
    else
        fprintf('Error: level = <"velocity", "acceleration">');
    end