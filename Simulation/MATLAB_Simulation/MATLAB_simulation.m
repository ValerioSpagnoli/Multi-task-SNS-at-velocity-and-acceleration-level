classdef MATLAB_simulation
    properties
        robot
        robot_model
        joint_positions
        directional_error
        T 
        points
    end

    methods
        function self = MATLAB_simulation(robot, joint_positions, directional_error, simulation_step, points)
            
            clear figure;

            self.robot = robot;
            if strcmp(robot.name, 'KUKA_LBR_IIWA_7_R800')
                self.robot_model = importrobot('iiwa7.urdf'); 
                self.robot_model.DataFormat = 'row';
            elseif strcmp(robot.name, 'KUKA_LBR_IV')
                disp('KUKA LBR IV does not exist in Matlab System Toolbox. The simuluation will be done showing the end effector position only.');
                self.robot_model = NaN;
            end

            self.joint_positions = joint_positions; 
            self.directional_error = directional_error;
            self.T = simulation_step;
            self.points = points;
            
            if ~isnan(joint_positions)
                self.run_3D_simulation();
            else
                disp('Cannot run 3D simulation because joint position is NaN.')
                return;
            end
            
            if ~isnan(directional_error)
                self.plot_directional_error();
            else
                disp('Cannot plot directional error because is NaN.')
                return;
            end
    
        end

        %% main function
        function run_3D_simulation(self)
            robot_arm = self.robot_model;
            ndof = self.robot.ndof;
            
            bounds_min_position = self.robot.bounds_position(1,:);
            bounds_max_position = self.robot.bounds_position(2,:);
                
            framerate = 15;
            r = rateControl(framerate);
            t_init = 1*self.T;
            t_final = 10000*self.T;
            num_frames = t_final*framerate;
            
            time = [0,linspace(t_init,t_final,size(self.joint_positions,2)-1)];
            qInterp = pchip(time,self.joint_positions,linspace(t_init,t_final,num_frames))';
            
            gripperPosition = zeros(num_frames,3);
            for k = 1:num_frames
                gripperPosition(k,:) = self.robot.get_ee_position(qInterp(k,:)');
            end
            

            figure;
            if isobject(robot_arm)
                show(robot_arm, self.joint_positions(:,1)', 'PreservePlot', false);
                hold on
            end
        
            points_x = zeros(1, length(self.points));
            points_y = zeros(1, length(self.points));
            points_z = zeros(1, length(self.points));
            for i=1:size(self.points,2)
                points_x(i) = self.points(1,i);
                points_y(i) = self.points(2,i);
                points_z(i) = self.points(3,i);
            end
        
            plot3(points_x, points_y, points_z, '-square', 'Color', 'r');
            hold on
            grid on
            p = plot3(gripperPosition(1,1), gripperPosition(1,2), gripperPosition(1,3), 'Color', 'black');
            hold on
            grid on
            
            pause(1);
            
            for k = 1:size(qInterp,1)
                if isobject(robot_arm)
                    show(robot_arm, qInterp(k,:), 'PreservePlot', false);
                end

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

        %% Plot directional error
        function plot_directional_error(self)
            t_init = 0.001;
            t_final = length(self.directional_error)/1000;
            time = [0,linspace(t_init,t_final,size(self.directional_error,2)-1)];

            figure;
            plot(time, self.directional_error);
            grid on
            hold on
            title('Directional error');
            xlabel('Time [s]');
            ylabel('Directional error [rad]');
        end


    end
end