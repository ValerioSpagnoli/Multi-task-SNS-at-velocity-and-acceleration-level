classdef MATLAB_simulation
    properties
        robot
        robot_model
        joint_positions
        directional_error
        x_dot_d
        x_ddot_d
        versors
        T 
        points
    end

    methods
        % function self = MATLAB_simulation(robot, joint_positions, directional_error, simulation_step, points)
        function self = MATLAB_simulation(sim)
            
            clear figure;

            self.robot = sim.robot;
            if strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800')
                self.robot_model = importrobot('iiwa7.urdf'); 
                self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'KUKA_LBR_IV')
                disp('KUKA LBR IV does not exist in Matlab System Toolbox. The simuluation will be done showing the end effector position only.');
                self.robot_model = importrobot('iiwa7.urdf');
                self.robot_model.DataFormat = 'row';
            end

            self.joint_positions = sim.joint_positions; 
            self.directional_error = sim.directional_error;
            self.x_dot_d = sim.x_dot_d;
            self.x_ddot_d = sim.x_ddot_d;
            self.versors = sim.versors;
            self.T = sim.simulation_step;
            self.points = sim.path;
            
            if ~isnan(sim.joint_positions)
                self.run_3D_simulation();
            else
                disp('Cannot run 3D simulation because joint position is NaN.')
                return;
            end
            
            if ~isnan(sim.directional_error)
                self.plot_directional_error();
            else
                disp('Cannot plot directional error because is NaN.')
                return;
            end

            if ~isnan(sim.x_dot_d)
                self.plot_x_dot_d();
            else
                disp('Cannot plot x_dot_d because is NaN.')
                return;
            end

            if ~isnan(sim.x_ddot_d)
                self.plot_x_ddot_d();
            else
                disp('Cannot plot x_ddot_d because is NaN.')
                return;
            end

%             points_x = zeros(1, length(self.points));
%             points_y = zeros(1, length(self.points));
%             points_z = zeros(1, length(self.points));
%             for i=1:size(self.points,2)
%                 points_x(i) = self.points(1,i);
%                 points_y(i) = self.points(2,i);
%                 points_z(i) = self.points(3,i);
%             end
% 
%             plot3(points_x, points_y, points_z, '-square', 'Color', 'r');
%             hold on
%             grid on
    
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
            for i=1:size(self.versors,2)
                v = self.versors{i};
                pi = v(:,1);
                pf = v(:,2); 
                quiver3(pi(1), pi(2), pi(3), pf(1), pf(2), pf(3));
            end


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

        function plot_x_dot_d(self)
            t_init = 0.001;
            t_final = length(self.x_dot_d)/1000;
            time = [0,linspace(t_init,t_final,size(self.x_dot_d,2)-1)];

            figure;
            plot(time, self.x_dot_d);
            grid on
            hold on
            title('x_dot_D');
            xlabel('Time [s]');
            ylabel('x_dot_d [m/s]');
        end

        function plot_x_ddot_d(self)
            t_init = 0.001;
            t_final = length(self.x_ddot_d)/1000;
            time = [0,linspace(t_init,t_final,size(self.x_ddot_d,2)-1)];

            figure;
            plot(time, self.x_ddot_d);
            grid on
            hold on
            title('x_ddot_d');
            xlabel('Time [s]');
            ylabel('x_ddot_d [m/s^2]');
            
            

        end

    end
end