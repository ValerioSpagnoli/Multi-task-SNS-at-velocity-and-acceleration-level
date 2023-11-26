classdef MATLAB_simulation
    properties
        robot
        robot_model
        bounds
        joints_positions
        joints_velocities
        joints_accelerations
        ee_positions
        ee_velocities
        link_1_positions
        link_1_velocities
        elbow_positions
        elbow_velocities
        directional_errors
        T 
        points
        folder_path
        title_font_size
        label_font_size
        legend_font_size
    end

    methods        
        function self = MATLAB_simulation(sim)
            
            clear figure;
            self.folder_path = '/Users/valeriospagnoli/Documents/MATLAB/Robotics/MultiTask_SNS/Simulation/Results/KUKA_LBR_IV/Segment/task1_task2/Plots/';
            self.title_font_size = 30;
            self.label_font_size = 22;
            self.legend_font_size = 18;

            self.robot = sim.robot;
            if strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800')
                self.robot_model = importrobot('iiwa7.urdf'); 
                self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'KUKA_LBR_IV')                
                %[self.robot_model] = importrobot('KUKA_LBR_IV.urdf');                             
                %self.robot_model.DataFormat = 'row';
            elseif strcmp(self.robot.name, 'Planar4R')
                self.robot_model = importrobot('Planar4R.urdf', 'urdf'); 
                self.robot_model.DataFormat = 'row';
            end

            self.bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
            self.joints_positions = sim.joints_positions; 
            self.joints_velocities = sim.joints_velocities; 
            self.joints_accelerations = sim.joints_accelerations; 
            self.ee_positions = sim.ee_positions;
            self.ee_velocities = sim.ee_velocities;
            if strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                self.elbow_positions = sim.elbow_positions;
                self.elbow_velocities = sim.elbow_velocities;
            elseif strcmp(self.robot.name, 'Planar4R')
                self.link_1_positions = sim.link_1_positions;
                self.link_1_velocities = sim.link_1_velocities;
            end

            self.directional_errors = sim.directional_errors;
            self.T = sim.simulation_step;
            self.points = sim.path;            
            
            if ~isnan(sim.joints_positions)
                self.run_3D_simulation();
            else
                disp('Cannot run 3D simulation because joints position vector is NaN.');
                return;
            end

            if ~isnan(self.joints_positions)
                self.plot_joints_positions();
            else
                disp('Cannot plot joints positions vector because is NaN.');
                return;
            end

            if ~isnan(self.joints_velocities)
                self.plot_joints_velocities();
            else
                disp('Cannot plot joints velocities vector because is NaN.');
                return;
            end

            if ~isnan(self.joints_accelerations)
                self.plot_joints_accelerations();
            else
                disp('Cannot plot joints accelerations vector because is NaN.');
                return;
            end
            
            if ~isnan(sim.directional_errors)
                self.plot_directional_error();
            else
                disp('Cannot plot directional error vector because is NaN.');
                return;
            end     

            if ~isnan(self.joints_positions)
                self.plot_link_position();
            else 
                disp('Cannot plot link cartisian positions because joints positions vector is NaN');
                return;
            end

            self.plot_link_velocities()

            self.plot_table();
        end

        %% main function
        function run_3D_simulation(self)
            robot_arm = self.robot_model;

            f = figure('Position', [100, 100, 1000, 800]);            

            framerate = 20;
            r = rateControl(framerate);
            t_init = 1*self.T;            
            t_final = double(length(self.joints_positions)*self.T);
            num_frames = int32(t_final*framerate)+1;       
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];

            interpolated_joints_positions = pchip(time,self.joints_positions,linspace(t_init,t_final,num_frames))';
            interpolated_ee_positions = pchip(time,self.ee_positions,linspace(t_init,t_final,num_frames))';

            if strcmp(self.robot.name, 'Planar4R')
                interpolated_link_1_positions = pchip(time,self.link_1_positions,linspace(t_init,t_final,num_frames))';                
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                interpolated_elbow_positions = pchip(time,self.elbow_positions,linspace(t_init,t_final,num_frames))';                
            end
            
            if isobject(robot_arm)
                show(robot_arm, self.joints_positions(:,1)', 'PreservePlot', false); hold on;
            end

            if strcmp(self.robot.name, 'Planar4R')
                if isobject(robot_arm)
                    az = 0; 
                    el = 90;  
                    
                    zoomFactor = 1;   
                    zoomPoint = [0.7; 0; 0];
                    xlim([zoomPoint(1) - zoomFactor, zoomPoint(1) + zoomFactor]);
                    ylim([zoomPoint(2) - zoomFactor, zoomPoint(2) + zoomFactor]);
                    zlim([0, zoomPoint(3) + zoomFactor]);
    
                    [x, y, z] = meshgrid(-5:0.2:5, -5:0.2:5, 0:0);
                    surf(x, y, z, 'FaceColor','0.1,0.1,0.1', 'FaceAlpha','0.1'); hold on; 
                    view([az, el]);
                else
                    az = 0; 
                    el = 90;  
                    
                    zoomFactor = 1;   
                    zoomPoint = [0.7; 0; 0];
                    xlim([zoomPoint(1) - zoomFactor, zoomPoint(1) + zoomFactor]);
                    ylim([zoomPoint(2) - zoomFactor, zoomPoint(2) + zoomFactor]);
                    zlim([0, zoomPoint(3) + zoomFactor]);
    
                    [x, y, z] = meshgrid(0:0.2:1.2, -0.6:0.2:0.6, 0:0);
                    surf(x, y, z, 'FaceColor','1,1,1', 'FaceAlpha','1'); hold on; 
                    view([az, el]);
                end
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                if isobject(robot_arm)
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
                    view([az, el]);
                else
                    az = 0; 
                    el = 45;  
                    
                    zoomFactor = 0.8;   
                    zoomPoint = [0.1; 0.3; 0.7];
                    xlim([zoomPoint(1) - zoomFactor, zoomPoint(1) + zoomFactor]);
                    ylim([zoomPoint(2) - zoomFactor, zoomPoint(2) + zoomFactor]);
                    zlim([0, zoomPoint(3) + zoomFactor]);
   
                    view([az, el]);
                end
            end
            

            points_x = zeros(1, length(self.points));
            points_y = zeros(1, length(self.points));
            points_z = zeros(1, length(self.points));
            for i=1:size(self.points,2)
                points_x(i) = self.points(1,i);
                points_y(i) = self.points(2,i);
                points_z(i) = self.points(3,i);
            end
            plot3(points_x, points_y, points_z, '-square', 'Color', 'r'); hold on;
                    
            p1 = plot3(self.ee_positions(1,1), self.ee_positions(1,2), self.ee_positions(1,3), 'Color', 'black', 'LineWidth',1.5); hold on;

            if strcmp(self.robot.name, 'Planar4R')
                p2 = plot3(interpolated_link_1_positions(1,1), interpolated_link_1_positions(1,2), interpolated_link_1_positions(1,3), 'Color', 'blue', 'LineWidth',1.5); hold on;                               
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                p2 = plot3(interpolated_elbow_positions(1,1), interpolated_elbow_positions(1,2), interpolated_elbow_positions(1,3), 'Color', 'blue', 'LineWidth',1.5); hold on;                
            end
            grid on;

            pause(2);
            
            for k = 1:size(interpolated_joints_positions,1)
                if isobject(robot_arm)
                    show(robot_arm, interpolated_joints_positions(k,:), 'PreservePlot', false);                    
                end

                bounds_min_position = self.bounds{1}(1,:);
                bounds_max_position = self.bounds{1}(2,:);
                
                for i=1:self.robot.n
                    if interpolated_joints_positions(k,i) > bounds_max_position(i) || interpolated_joints_positions(k,i) < bounds_min_position(i)
                        fprintf('Joint %d out of bounds\n', i);
                    end
                end
              
                p1.XData(k) = interpolated_ee_positions(k,1);
                p1.YData(k) = interpolated_ee_positions(k,2);
                p1.ZData(k) = interpolated_ee_positions(k,3);

                if strcmp(self.robot.name, 'Planar4R')
                    p2.XData(k) = interpolated_link_1_positions(k,1);
                    p2.YData(k) = interpolated_link_1_positions(k,2);
                    p2.ZData(k) = interpolated_link_1_positions(k,3);
                                       
                elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                    p2.XData(k) = interpolated_elbow_positions(k,1);
                    p2.YData(k) = interpolated_elbow_positions(k,2);
                    p2.ZData(k) = interpolated_elbow_positions(k,3);                    
                end
                waitfor(r);
            end
            
            if strcmp(self.robot.name, 'Planar4R')
                title('End-effector and link_1 trajectories', 'FontSize',self.title_font_size);
                xlabel('X [m]', 'FontSize',self.label_font_size);
                ylabel('Y [m]', 'FontSize',self.label_font_size);
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                title('End-effector and elbow trajectories', 'FontSize',26);
                xlabel('X [m]', 'FontSize',16);
                ylabel('Y [m]', 'FontSize',16);
                zlabel('Z [m]', 'FontSize',16);
            end

            if isobject(robot_arm)
                exportgraphics(f, strcat(self.folder_path,'path_robot.png'));
            else
                exportgraphics(f, strcat(self.folder_path,'path.png'));
            end
            hold off
        end

        %% Plot directional error
        function plot_directional_error(self)
            f = figure('Position', [100, 100, 1000, 300]);            

            t_init = self.T;
            t_final = length(self.directional_errors)*self.T;
            time = [0,linspace(t_init,t_final,size(self.directional_errors,2)-1)];
            
            plot(time, self.directional_errors);
            grid on
            hold on
            title('Directional errors', 'FontSize',self.title_font_size);
            xlabel('Time [s]', 'FontSize',self.label_font_size);
            ylabel('Directional errors [rad]', 'FontSize',self.label_font_size);            
            exportgraphics(f, strcat(self.folder_path,'directional_errors.png'));
        end
        
        %% Plot joints positions
        function plot_joints_positions(self)               
            f = figure('Position', [100, 100, 1000, 300]);

            t_init = self.T;
            t_final = length(self.joints_positions)*self.T;        
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];            

            for i=1:self.robot.n
                plot(time, self.joints_positions(i,:), 'DisplayName', sprintf('Joint_%d position',i));hold on;                
            end
            
            grid on
            title('Joints positions', 'FontSize',self.title_font_size);
            xlabel('Time [s]', 'FontSize',self.label_font_size);
            ylabel('Joints positions [rad]', 'FontSize',self.label_font_size);
            legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);            
            exportgraphics(f, strcat(self.folder_path,'joints_positions.png'));
        end
        
        %% Plot joints velocities
        function plot_joints_velocities(self)   
            f = figure('Position', [100, 100, 1000, 300]);         
            
            t_init = self.T;
            t_final = length(self.joints_velocities)*self.T;
            time = [0,linspace(t_init,t_final,size(self.joints_velocities,2)-1)];
            
            for i=1:self.robot.n
                plot(time, self.joints_velocities(i,:), 'DisplayName', sprintf('Joint_%d velocity',i));hold on;
            end
            
            grid on
            title('Joints velocities', 'FontSize',self.title_font_size);
            xlabel('Time [s]', 'FontSize',self.label_font_size);
            ylabel('Joints velocities [rad/s]', 'FontSize',self.label_font_size);
            legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);            
            exportgraphics(f, strcat(self.folder_path,'joints_velocities.png'));
        end

        %% Plot joints accelerations
        function plot_joints_accelerations(self) 
            f = figure('Position', [100, 100, 1000, 300]);
            
            t_init = self.T;
            t_final = length(self.joints_accelerations)*self.T;
            time = [0,linspace(t_init,t_final,size(self.joints_accelerations,2)-1)];

            for i=1:self.robot.n                
                plot(time, self.joints_accelerations(i,:), 'DisplayName', sprintf('Joint_%d acceleration',i));hold on;
            end
            
            grid on
            title('Joints accelerations', 'FontSize',self.title_font_size);
            xlabel('Time [s]', 'FontSize',self.label_font_size);
            ylabel('Joints accelerations [rad/s^2]', 'FontSize',self.label_font_size);
            legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);            
            exportgraphics(f, strcat(self.folder_path,'joints_accelerations.png'));
        end


        %% Plot elbow position - link 1 position
        function plot_link_position(self) 
            f = figure('Position', [100, 100, 1000, 300]);
            
            t_init = self.T;
            t_final = length(self.joints_positions)*self.T;
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];
                                    
            if strcmp(self.robot.name, 'Planar4R')               
                plot(time, self.link_1_positions(1,:), 'DisplayName', 'Link 1 position x');hold on;                
                plot(time, self.link_1_positions(2,:), 'DisplayName', 'Link 1 position y');hold on;                
                plot(time, self.link_1_positions(3,:), 'DisplayName', 'Link 1 position z');hold on;                                
                
                grid on
                title('Link 1 positions', 'FontSize',self.title_font_size);
                xlabel('Time [s]', 'FontSize',self.label_font_size);
                ylabel('Link 1 positions [m]', 'FontSize',self.label_font_size);
                legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);                
                exportgraphics(f, strcat(self.folder_path,'link1_positions.png'));

            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                plot(time, self.elbow_positions(1,:), 'DisplayName', 'Elbow position x');hold on;                
                plot(time, self.elbow_positions(2,:), 'DisplayName', 'Elbow position y');hold on;                
                plot(time, self.elbow_positions(3,:), 'DisplayName', 'Elbow position z');hold on;                
                           
                grid on
                title('Elbow positions', 'FontSize',self.title_font_size);
                xlabel('Time [s]', 'FontSize',self.label_font_size);
                ylabel('Elbow positions [m]', 'FontSize',self.label_font_size);
                legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);                
                exportgraphics(f, strcat(self.folder_path,'elbow_positions.png'));
            end            
        end

         %% Plot elbow position - link 1 position
        function plot_link_velocities(self) 
            f = figure('Position', [100, 100, 1000, 300]);              
            
            t_init = self.T;
            t_final = length(self.joints_positions)*self.T;
            time = [0,linspace(t_init,t_final,size(self.joints_positions,2)-1)];           
            
            if strcmp(self.robot.name, 'Planar4R')               
                plot(time, self.link_1_velocities(1,:), 'DisplayName', 'Link 1 velocity x');hold on;                
                plot(time, self.link_1_velocities(2,:), 'DisplayName', 'Link 1 velocity y');hold on;                
                plot(time, self.link_1_velocities(3,:), 'DisplayName', 'Link 1 velocity z');hold on;                                
                
                grid on
                title('Link 1 velocities', 'FontSize',self.title_font_size);
                xlabel('Time [s]', 'FontSize',self.label_font_size);
                ylabel('Link 1 velocities [m/s]', 'FontSize',self.label_font_size);
                legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);                
                exportgraphics(f, strcat(self.folder_path,'link1_velocities.png'));

            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                plot(time, self.elbow_velocities(1,:), 'DisplayName', 'Elbow velocity x');hold on;                
                plot(time, self.elbow_velocities(2,:), 'DisplayName', 'Elbow velocity y');hold on;                
                plot(time, self.elbow_velocities(3,:), 'DisplayName', 'Elbow velocity z');hold on;                
                           
                grid on
                title('Elbow velocities', 'FontSize',self.title_font_size);
                xlabel('Time [s]', 'FontSize',self.label_font_size);
                ylabel('Elbow velocities [m/s]', 'FontSize',self.label_font_size);
                legend('Location', 'best', 'Orientation', 'vertical', 'FontSize',self.legend_font_size);                
                exportgraphics(f, strcat(self.folder_path,'elbow_velocities.png'));
            end
        end


        %% Plot table
        %  The table contain the following metrics:
        %  - total time
        %  - average directional error
        %  - average elbow velocity (for KUKAs) or 
        %    average link_1 velocities (for Planar4R)
        %  - average of the absolute value of the elbow position (for KUKAs) or 
        %    average of the absolute value of the link_1 (for Planar4R)
            
        function plot_table(self)                            
            % Total time
            Total_time = length(self.joints_accelerations)*self.T;
            
            % average direction error
            average_directional_error = mean(self.directional_errors);

            % average joint position and velocities
            if strcmp(self.robot.name, 'Planar4R')
                avg_norm_link_1_velocities = 0;
                avg_link_1_x_positions = 0;
                avg_link_1_y_positions = 0;
                for i=1:length(self.link_1_velocities)
                    avg_norm_link_1_velocities = avg_norm_link_1_velocities + norm(self.link_1_velocities(:,i));
                    avg_link_1_x_positions = avg_link_1_x_positions + self.link_1_positions(1,i);
                    avg_link_1_y_positions = avg_link_1_y_positions + self.link_1_positions(2,i);
                end
                avg_norm_link_1_velocities = avg_norm_link_1_velocities/length(self.link_1_velocities);
                avg_link_1_x_positions = avg_link_1_x_positions/length(self.link_1_positions);   
                avg_link_1_y_positions = avg_link_1_y_positions/length(self.link_1_positions);   
                
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
                avg_norm_elbow_velocities = 0;
                avg_elbow_x_positions = 0;
                avg_elbow_y_positions = 0;
                avg_elbow_z_positions = 0;
                for i=1:length(self.elbow_velocities)
                    avg_norm_elbow_velocities = avg_norm_elbow_velocities + norm(self.elbow_velocities(:,i));
                    avg_elbow_x_positions = avg_elbow_x_positions + self.elbow_positions(1,i);
                    avg_elbow_y_positions = avg_elbow_y_positions + self.elbow_positions(2,i);
                    avg_elbow_z_positions = avg_elbow_z_positions + self.elbow_positions(3,i);
                end
                avg_norm_elbow_velocities = avg_norm_elbow_velocities/length(self.elbow_velocities);
                avg_elbow_x_positions = avg_elbow_x_positions/length(self.elbow_positions);   
                avg_elbow_y_positions = avg_elbow_y_positions/length(self.elbow_positions); 
                avg_elbow_z_positions = avg_elbow_z_positions/length(self.elbow_positions); 
            end
            
            fprintf('------------- TABLE -------------\n');
            fprintf('Total time [s]:                    ');disp(Total_time);
            fprintf('Average directional error: [rad]   ');disp(average_directional_error);
            if strcmp(self.robot.name, 'Planar4R')
            fprintf('Average link 1 norm velocities: [m/s]   ');disp(avg_norm_link_1_velocities);
            fprintf('Average link 1 x-position: [m]     ');disp(avg_link_1_x_positions);
            fprintf('Average link 1 y-position: [m]     ');disp(avg_link_1_y_positions);
            elseif strcmp(self.robot.name, 'KUKA_LBR_IIWA_7_R800') || strcmp(self.robot.name, 'KUKA_LBR_IV')
            fprintf('Average elbow norm velocities:  [m/s]   ');disp(avg_norm_elbow_velocities);
            fprintf('Average elbow x-position:  [m]     ');disp(avg_elbow_x_positions);
            fprintf('Average elbow y-position:  [m]     ');disp(avg_elbow_y_positions);
            fprintf('Average elbow z-position:  [m]     ');disp(avg_elbow_z_positions);
            end
        end

    end
end