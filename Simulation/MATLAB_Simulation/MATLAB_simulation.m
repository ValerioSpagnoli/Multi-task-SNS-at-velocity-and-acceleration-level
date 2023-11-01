classdef MATLAB_simulation
    properties
        robot
        robot_model
        bounds
        joints_positions
        joints_velocities
        joints_accelerations
        directional_errors
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
                [self.robot_model] = importrobot('KUKA_LBR_IV.urdf');                             
                self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'Planar4R')
                self.robot_model = importrobot('Planar4R.urdf', 'urdf'); 
                self.robot_model.DataFormat = 'row';
            end

            self.bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
            self.joints_positions = sim.joints_positions; 
            self.joints_velocities = sim.joints_velocities; 
            self.joints_accelerations = sim.joints_accelerations; 
            self.directional_errors = sim.directional_errors;
            self.T = sim.simulation_step;
            self.points = sim.path;
            
            if ~isnan(sim.joints_positions)
                self.run_3D_simulation();
            else
                disp('Cannot run 3D simulation because joints position vector is NaN.')
                return;
            end

            if ~isnan(self.joints_positions)
                self.plot_joints_positions();
            else
                disp('Cannot plot joints positions vector because is NaN.')
                return;
            end

            if ~isnan(self.joints_velocities)
                self.plot_joints_velocities();
            else
                disp('Cannot plot joints velocities vector because is NaN.')
                return;
            end

            if ~isnan(self.joints_accelerations)
                self.plot_joints_accelerations();
            else
                disp('Cannot plot joints accelerations vector because is NaN.')
                return;
            end
            
            if ~isnan(sim.directional_errors)
                self.plot_directional_error();
            else
                disp('Cannot plot directional error vector because is NaN.')
                return;
            end     
        end

        %% main function
        function run_3D_simulation(self)
            robot_arm = self.robot_model;
            n = self.robot.n;
            
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


            if strcmp(self.robot.name, 'Planar4R')
                link_1_positions = zeros(num_frames,3);
                link_2_positions = zeros(num_frames,3);
                link_3_positions = zeros(num_frames,3);       

                for k = 1:num_frames
                   link_1_positions(k,:) = self.robot.get_link_1_position(qInterp(k,:)');
                   link_2_positions(k,:) = self.robot.get_link_2_position(qInterp(k,:)');
                   link_3_positions(k,:) = self.robot.get_link_3_position(qInterp(k,:)');
                end

            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                elbow_position = zeros(num_frames,3);
                for k = 1:num_frames
                    elbow_position(k,:) = self.robot.get_elbow_position(qInterp(k,:)');
                end
            end

            
            figure('Position', [100, 100, 1000, 800]);

            if isobject(robot_arm)
                show(robot_arm, self.joints_positions(:,1)', 'PreservePlot', false); hold on;
            end

            if strcmp(self.robot.name, 'Planar4R')
                az = 0; 
                el = 90;  
                
                zoomFactor = 1;   
                zoomPoint = [0.7; 0; 0];
                xlim([zoomPoint(1) - zoomFactor, zoomPoint(1) + zoomFactor]);
                ylim([zoomPoint(2) - zoomFactor, zoomPoint(2) + zoomFactor]);
                zlim([0, zoomPoint(3) + zoomFactor]);

                [x, y, z] = meshgrid(-5:0.2:5, -5:0.2:5, 0:0);
                surf(x, y, z, 'FaceColor','0.1,0.1,0.1', 'FaceAlpha','0.1'); hold on; 
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                az = 80;
                el = 15;
                
                xlim([-1 1]);
                ylim([-1 1]);
                zlim([0 1.5]);

                zoomFactor = 0.8;   
                zoomPoint = [0.1; 0.3; 0.7];
                xlim([zoomPoint(1) - zoomFactor, zoomPoint(1) + zoomFactor]);
                ylim([zoomPoint(2) - zoomFactor, zoomPoint(2) + zoomFactor]);
                zlim([0, zoomPoint(3) + zoomFactor]);
                [x, y, z] = meshgrid(-5:0.2:5, -5:0.2:5, 0:0);
                surf(x, y, z, 'FaceColor','0.2,0.2,0.2', 'FaceAlpha','0.2'); hold on; 
            end
            view([az, el]);

            points_x = zeros(1, length(self.points));
            points_y = zeros(1, length(self.points));
            points_z = zeros(1, length(self.points));
            for i=1:size(self.points,2)
                points_x(i) = self.points(1,i);
                points_y(i) = self.points(2,i);
                points_z(i) = self.points(3,i);
            end
            plot3(points_x, points_y, points_z, '-square', 'Color', 'r'); hold on;
        
            p1 = plot3(end_effector_position(1,1), end_effector_position(1,2), end_effector_position(1,3), 'Color', 'black', 'LineWidth',1.5); hold on;

            if strcmp(self.robot.name, 'Planar4R')
                p2 = plot3(link_1_positions(1,1), link_1_positions(1,2), link_1_positions(1,3), 'Color', 'blue', 'LineWidth',1.5); hold on;               
                p3 = plot3(link_2_positions(1,1), link_2_positions(1,2), link_2_positions(1,3), 'Color', 'red'); hold on;
                p4 = plot3(link_3_positions(1,1), link_3_positions(1,2), link_3_positions(1,3), 'Color', 'magenta'); hold on;
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                p2 = plot3(elbow_position(1,1), elbow_position(1,2), elbow_position(1,3), 'Color', 'blue', 'LineWidth',1.5); hold on;
            end
            grid on;

            pause(10);
            
            for k = 1:size(qInterp,1)
                if isobject(robot_arm)
                    show(robot_arm, qInterp(k,:), 'PreservePlot', false);
                end

                q_k = qInterp(k,:);

                for i=1:n
                    if q_k(i) > bounds_max_position(i)
                        fprintf('Joint %d over bounds_max_position: q_k(%d) = %f, bounds_max_position(%d) = %f\n', i, i, q_k(i), i, bounds_max_position(i));
                    elseif q_k(i) < bounds_min_position(i)
                        fprintf('Joint %d under bounds_min_position: q_k(%d) = %f, bounds_min_position(%d) = %f\n', i, i, q_k(i), i, bounds_min_position(i));
                    end
                end

                p1.XData(k) = end_effector_position(k,1);
                p1.YData(k) = end_effector_position(k,2);
                p1.ZData(k) = end_effector_position(k,3);


                if strcmp(self.robot.name, 'Planar4R')
                    p2.XData(k) = link_1_positions(k,1);
                    p2.YData(k) = link_1_positions(k,2);
                    p2.ZData(k) = link_1_positions(k,3);

                    p3.XData(k) = link_2_positions(k,1);
                    p3.YData(k) = link_2_positions(k,2);
                    p3.ZData(k) = link_2_positions(k,3);
                    
                    p4.XData(k) = link_3_positions(k,1);
                    p4.YData(k) = link_3_positions(k,2);
                    p4.ZData(k) = link_3_positions(k,3);
                    
                elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                    p2.XData(k) = elbow_position(k,1);
                    p2.YData(k) = elbow_position(k,2);
                    p2.ZData(k) = elbow_position(k,3);
                end
                waitfor(r);
            end
            hold off
        end

        %% Plot directional error
        function plot_directional_error(self)
            t_init = 0.001;
            t_final = length(self.directional_errors)/1000;
            time = [0,linspace(t_init,t_final,size(self.directional_errors,2)-1)];

            figure('Position', [100, 100, 1000, 800]);
            plot(time, self.directional_errors);
            grid on
            hold on
            title('Directional errors');
            xlabel('Time [s]');
            ylabel('Directional errors [rad]');
        end
        
        %% Plot joints positions
        function plot_joints_positions(self)               
            joints_pos = cell(1,self.robot.n);
            joints_pos{1} = [];
            joints_pos{2} = [];
            joints_pos{3} = [];
            joints_pos{4} = [];           
            for i=1:self.robot.n
                for j=1:length(self.joints_positions)
                    joints_pos{i} = [joints_pos{i}, self.joints_positions(i,j)];
                end
            end

            t_init = 0.001;
            t_final = length(self.joints_positions)/1000;        
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];            

            figure('Position', [100, 100, 1000, 800]);
            for i=1:self.robot.n
                plot(time, joints_pos{i}, 'DisplayName', sprintf('Joint_%d position',i));hold on;                
            end
            
            grid on
            title('Joints positions');
            xlabel('Time [s]');
            ylabel('Joints positions [rad]');

            legend('Location', 'best', 'Orientation', 'vertical');
        end
        
        %% Plot joints velocities
        function plot_joints_velocities(self)   
            joints_vel = cell(1,self.robot.n);
            joints_vel{1} = [];
            joints_vel{2} = [];
            joints_vel{3} = [];
            joints_vel{4} = [];           
            for i=1:self.robot.n
                for j=1:length(self.joints_velocities)
                    joints_vel{i} = [joints_vel{i}, self.joints_velocities(i,j)];
                end
            end

            t_init = 0.001;
            t_final = length(self.joints_velocities)/1000;
            time = [0,linspace(t_init,t_final,size(self.joints_velocities,2)-1)];

            figure('Position', [100, 100, 1000, 800]);
            for i=1:self.robot.n
                plot(time, joints_vel{i}, 'DisplayName', sprintf('Joint_%d velocity',i));hold on;
            end
            
            grid on
            title('Joints velocities');
            xlabel('Time [s]');
            ylabel('Joints velocities [rad/s]');

            legend('Location', 'best', 'Orientation', 'vertical');
        end

        %% Plot joints accelerations
        function plot_joints_accelerations(self) 
            joints_acc = cell(1,self.robot.n);
            joints_acc{1} = [];
            joints_acc{2} = [];
            joints_acc{3} = [];
            joints_acc{4} = [];           
            for i=1:self.robot.n
                for j=1:length(self.joints_accelerations)                    
                    joints_acc{i} = [joints_acc{i}, self.joints_accelerations(i,j)];
                end
            end

            t_init = 0.001;
            t_final = length(self.joints_accelerations)/1000;
            time = [0,linspace(t_init,t_final,size(self.joints_accelerations,2)-1)];

            figure('Position', [100, 100, 1000, 800]);
            for i=1:self.robot.n                
                plot(time, joints_acc{i}, 'DisplayName', sprintf('Joint_%d acceleration',i));hold on;
            end
            
            grid on
            title('Joints accelerations');
            xlabel('Time [s]');
            ylabel('Joints accelerations [rad/s^2]');

            legend('Location', 'best', 'Orientation', 'vertical');
        end
    end
end