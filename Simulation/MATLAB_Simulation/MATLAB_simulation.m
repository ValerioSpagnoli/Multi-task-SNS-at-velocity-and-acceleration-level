classdef MATLAB_simulation
    properties
        robot
        robot_model
        joints_positions
        directional_errors
        elbow_positions
        elbow_velocities
        T 
        points
    end

    methods        
        function self = MATLAB_simulation(sim)
            
            clear figure;

            self.robot = sim.robot;
            if strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800')
                self.robot_model = importrobot('iiwa7.urdf'); 
                self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'KUKA_LBR_IV')
                disp('KUKA LBR IV does not exist in Matlab System Toolbox. The simuluation will be done showing the end effector position only.');
                %self.robot_model = importrobot('iiwa7.urdf');
                %self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'Planar4R')
                self.robot_model = importrobot('Planar4R.urdf', 'urdf'); 
                self.robot_model.DataFormat = 'row';
            end

            self.joints_positions = sim.joints_positions; 
            self.directional_errors = sim.directional_errors;
            self.elbow_positions = sim.elbow_positions; 
            self.elbow_velocities = sim.elbow_velocities;
            self.T = sim.simulation_step;
            self.points = sim.path;
            
            if ~isnan(sim.joints_positions)
                self.run_3D_simulation();
            else
                disp('Cannot run 3D simulation because joints position vector is NaN.')
                return;
            end
            
            if ~isnan(sim.directional_errors)
                self.plot_directional_error();
            else
                disp('Cannot plot directional error vector because is NaN.')
                return;
            end

            if ~isnan(sim.elbow_velocities)
                self.plot_elbow_velocities();
            else
                disp('Cannot plot elbow velocities vector because is NaN.')
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
%     
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
            
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];
            qInterp = pchip(time,self.joints_positions,linspace(t_init,t_final,num_frames))';
            
            end_effector_position = zeros(num_frames,3);
            for k = 1:num_frames
                end_effector_position(k,:) = self.robot.get_ee_position(qInterp(k,:)');
            end

            elbow_position = zeros(num_frames,3);
%             for k = 1:num_frames
%                 elbow_position(k,:) = self.robot.get_elbow_position(qInterp(k,:)');
%             end

            
            figure;
            if isobject(robot_arm)
                show(robot_arm, self.joints_positions(:,1)', 'PreservePlot', false);
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
        
            p1 = plot3(end_effector_position(1,1), end_effector_position(1,2), end_effector_position(1,3), 'Color', 'black');            
            hold on
            p2 = plot3(elbow_position(1,1), elbow_position(1,2), elbow_position(1,3), 'Color', 'blue');
            hold on
            grid on
            
            pause(10);
            
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

                p1.XData(k) = end_effector_position(k,1);
                p1.YData(k) = end_effector_position(k,2);
                p1.ZData(k) = end_effector_position(k,3);

                p2.XData(k) = elbow_position(k,1);
                p2.YData(k) = elbow_position(k,2);
                p2.ZData(k) = elbow_position(k,3);

                waitfor(r);
            end
            hold off
        end

        %% Plot directional error
        function plot_directional_error(self)
            t_init = 0.001;
            t_final = length(self.directional_errors)/1000;
            time = [0,linspace(t_init,t_final,size(self.directional_errors,2)-1)];

            figure;
            plot(time, self.directional_errors);
            grid on
            hold on
            title('Directional errors');
            xlabel('Time [s]');
            ylabel('Directional errors [rad]');
        end
        
        %% Plot elbow velocies
        function plot_elbow_velocities(self)
            t_init = 0.001;
            t_final = length(self.elbow_velocities)/1000;
            time = [0,linspace(t_init,t_final,size(self.elbow_velocities,2)-1)];

            figure;
            plot(time, self.elbow_velocities);
            grid on
            hold on
            title('Elbow velocities');
            xlabel('Time [s]');
            ylabel('Elbow velocities [rad/s]');
        end
    end
end