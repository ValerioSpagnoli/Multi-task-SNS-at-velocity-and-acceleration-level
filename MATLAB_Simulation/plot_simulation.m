function plot_simulation(T, qs, robot, points)
    ndof = robot.ndof;
    robot_arm = robot.robot;

    bounds_min_position = robot.bounds_position(1,:);
    bounds_max_position = robot.bounds_position(2,:);

    bounds_min_velocity = robot.bounds_velocity(1,:);
    bounds_max_velocity = robot.bounds_velocity(2,:);

    bounds_min_acceleration = robot.bounds_acceleration(1,:);
    bounds_max_acceleration = robot.bounds_acceleration(2,:);

    q_0 = robot.q_0;

%     robot_arm = robot{1};
%     bounds = robot{2};
%     ndof = robot{3};
%     q_0 = robot{4};
% 
%     bounds_min_position = bounds{1}(1,:);
%     bounds_max_position = bounds{1}(2,:);
% 
%     bounds_min_velocity = bounds{2}(1,:);
%     bounds_max_velocity = bounds{2}(2,:);
% 
%     bounds_min_acceleration = bounds{3}(1,:);
%     bounds_max_acceleration = bounds{3}(2,:);


    framerate = 15;
    r = rateControl(framerate);
    t_init = 1*T;
    t_final = 10000*T;
    num_frames = t_final*framerate;
    
    time_qs = [0,linspace(t_init,t_final,size(qs,2)-1)];
    qInterp = pchip(time_qs,qs,linspace(t_init,t_final,num_frames))';
    
    gripperPosition = zeros(num_frames,3);
    for k = 1:num_frames
        %gripperPosition(k,:) = tform2trvec(getTransform(robot_arm,qInterp(k,:),'iiwa_link_ee'));
        gripperPosition(k,:) = robot.get_ee_position(qInterp(k,:)');
    end
    
    figure;
    show(robot_arm, qs(:,1)', 'PreservePlot', false);
    hold on

    points_x = zeros(1, length(points));
    points_y = zeros(1, length(points));
    points_z = zeros(1, length(points));
    for i=1:length(points)
        points_x(i) = points(1,i);
        points_y(i) = points(2,i);
        points_z(i) = points(3,i);
    end
    plot3(points_x, points_y, points_z, '-square', 'Color', 'r')
    hold on

    pause(5);
    
    p = plot3(gripperPosition(1,1), gripperPosition(1,2), gripperPosition(1,3), 'Color', 'black');
    hold on

    for k = 1:size(qInterp,1)
        show(robot_arm, qInterp(k,:), 'PreservePlot', false);
        q_k = qInterp(k,:);
        for i=1:ndof
            if q_k(i) > bounds_max_position(i)
                fprintf('Joint %d over bounds_max_position: q_k(%d) = %f, bounds_max_position(%d) = %f\n', i, i, q_k(i), i, bounds_max_position(i));
            elseif q_k(i) < bounds_min_position(i)
                fprintf('Joint %d under bounds_min_position: q_k(%d) = %f, bounds_min_position(%d) = %f\n', i, i, q_k(i), i, bounds_min_position(i));
            end
        end
        p.XData(k) = gripperPosition(k,1);
        p.YData(k) = gripperPosition(k,2);
        p.ZData(k) = gripperPosition(k,3);
        waitfor(r);
    end
    hold off
end