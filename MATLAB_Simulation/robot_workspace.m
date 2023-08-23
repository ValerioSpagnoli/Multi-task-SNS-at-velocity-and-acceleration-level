function [px,py,pz] = robot_workspace(robot, bounds_position)

    bounds_min_position = bounds_position(1,:);
    bounds_max_position = bounds_position(2,:);

    i=1;
    px=zeros(1,189000);
    py=zeros(1,189000);
    pz=zeros(1,189000);
    
    for q1 = bounds_min_position(1):bounds_max_position(1)
    for q2 = bounds_min_position(2):bounds_max_position(2)
    for q3 = bounds_min_position(3):bounds_max_position(3)
    for q4 = bounds_min_position(4):bounds_max_position(4)
    for q5 = bounds_min_position(5):bounds_max_position(5)
    for q6 = bounds_min_position(6):bounds_max_position(6)
    for q7 = bounds_min_position(7):bounds_max_position(7)
        disp(i)
        T07 = getTransform(robot, [q1,q2,q3,q4,q5,q6,q7], 'iiwa_link_7');
        p = T07(1:3,4);
        px(i) = p(1);
        py(i) = p(2);
        pz(i) = p(3);
        i=i+1;
    end
    end
    end
    end
    end
    end
    end

    plot3(px, py, pz, '.', 'Color', 'c')
    grid on
    hold on