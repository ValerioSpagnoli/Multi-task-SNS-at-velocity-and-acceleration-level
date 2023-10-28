classdef Simulation_KUKA
    properties
        level
        robot
        robot_name
        q_0
        q_dot_0
        q_ddot_0
        simulation_step
        path
        epsilon
        joints_positions
        joints_velocities
        joints_accelerations
        directional_errors
    end

    methods
        %% Constructur
        function self = Simulation_KUKA(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon)

            format short

            % check parameters
            if isnan(level)
                self.level = 'velocity';
            else
                self.level = level;
            end

            if isnan(robot_name)
                self.robot_name = 'KUKA_LBR_IIWA_7_R800';
            else
                self.robot_name = robot_name;
            end

            if isnan(q_0)
                if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                    self.q_0 = [0; pi/4; -pi/4; pi/4; 0; 0; 0];
                elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                    self.q_0 = [0; pi/4; pi/4; pi/4; 0; 0; 0];
                end
            else
                self.q_0 = q_0;
            end

            if isnan(q_dot_0)                
                if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                    self.q_dot_0 = zeros(7,1);
                elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                    self.q_dot_0 = zeros(7,1);
                end
            else
                self.q_dot_0 = q_dot_0;    
            end

            if isnan(q_ddot_0)
                if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                    self.q_ddot_0 = zeros(7,1);
                elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                    self.q_ddot_0 = zeros(7,1);
                end
            else
                self.q_ddot_0 = q_ddot_0;
            end

            if isnan(simulation_step)
                self.simulation_step = 0.001;
            else
                self.simulation_step = simulation_step;
            end

            if isnan(epsilon)
                self.epsilon = 0.005;
            else
                self.epsilon = epsilon;
            end

            fprintf('Loading robot ... ');
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                self.robot = KUKA_LBR_IIWA7(self.q_0, self.q_dot_0, self.q_ddot_0);
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                self.robot = KUKA_LBR_IV(self.q_0, self.q_dot_0, self.q_ddot_0);
            end
            fprintf('done! Robot: %s\n', robot_name);
                     
            
            % compute path
            fprintf('Creating path ... ')
            self.path = self.create_hexagonal_path(2);
            % self.path = self.create_circular_path(1);
            fprintf('done! \n');

           
            % compute simulation
            fprintf('Start simulation ... \n\n')
            pause(1);
            if strcmp(level, 'velocity')
                [self.joints_positions, self.joints_velocities, self.directional_errors] = self.run_simulation_velocity_level();
            elseif strcmp(level, 'acceleration')
                [self.joints_positions, self.joints_velocities, self.joints_accelerations, self.directional_errors] = self.run_simulation_acceleration_level();
            end
            fprintf('Simulation ended. \n')
        end

        %% create path

        % Create hexagonal path
        function path = create_hexagonal_path(self, n_cycle)
            poly_0 = [[0; 0; 0.2], [0; 0.1732; 0.1], [0; 0.1732; -0.1], [0; 0; -0.2], [0; -0.1732; -0.1], [0; -0.1732; 0.1]];            
            
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                center = [0.1; 0.35; 0.8235];
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                center = [0.1; 0.35; 0.6235];
            end
                        
            points = [];
            for n=1:n_cycle
                for i=1:6
                    if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                        points = [points, [poly_0(:,i)]];
                    elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                        points = [points, [poly_0(:,i)]];
                    end
                end
            end
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                path = [points, poly_0(:,1)]+center;
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                path = [points, poly_0(:,1)]+center;
            end            
        end

        % Create circular path
        function path = create_circular_path(self, n_cycle)
            t = (0:0.2:10)'; 
            
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                center = [0.1; 0.35; 0.8235];
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                center = [0.1; 0.35; 0.6235];
            end
            radius = 0.2;
            theta = t*(2*pi/t(end));
            path = center + radius*[zeros(size(theta)) cos(theta) sin(theta)]';
        end

        %% run simulation velocity level
        function [joints_positions, joints_velocities, directional_errors] = run_simulation_velocity_level(self)

            n = self.robot.n;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
                            
            T = self.simulation_step;
            
            % Current configuration (time h)
            q_h = self.robot.q_0;
            q_dot_h = self.robot.q_dot_0;
            ee_position_h = self.robot.ee_position_0;            
            elbow_position_h = self.robot.elbow_position_0;

            % Jacobians of task 1, ee_position (time h)
            J1_h = self.robot.get_J_ee(q_h);
            J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);

            % Jacobians of task 2 (time h)
            J2_h = eye(n);
            J2_dot_h = zeros(7,7);

            % Jacobians of task 3, elbow_position (time h)
            J3_h = self.robot.get_J_eb(q_h);
            J3_dot_h = self.robot.get_J_dot_eb(q_h,q_dot_h);

            % Previous configuration (time h-1)
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;

            % Jacobian of task 1, ee_position (time h-1)
            J1_hm1 = self.robot.get_J_ee(q_hm1);                 
            J3_hm1 = self.robot.get_J_eb(q_hm1);   

            % Task 3:
            x3_dot_d_hm1 = zeros(3,1);      

            % Set PD controller gains of task 1
            Kp_1 = 2;
            Kd_1 = 0.05;
            
            % Set PD controller gains of task 3
            Kp_3 = 2;
            Kd_3 = 0.05;
       
            % All configurations during task
            joints_positions = [q_h];            
            joints_velocities = [q_dot_h];            
            directional_errors = [];
            
            k = 1;
            x_d = self.path(1:3,k);
            
            saturation_counter = 0;

            while true                
                if norm(ee_position_h - x_d) < 0.005
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                V_h = Kp_1*norm(x_d - ee_position_h) - Kd_1*norm(J1_hm1*q_dot_hm1);
                x1_dot_d_h = V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h));
                x1_ddot_d_h =  (x1_dot_d_h-J1_hm1*q_dot_hm1) / T;
                m1 = length(x1_ddot_d_h);               

                % TASK 2: self motion dumping
                q_dot_cs = -1000*q_h;
                m2 = length(q_dot_cs);

                % TASK 3: elbow
                elbow_target_point = [elbow_position_h(1); 0; elbow_position_h(3)];                
                if(norm(elbow_target_point-elbow_position_h)==0)
                   x3_dot_d_h = zeros(3,1);
                   x3_ddot_d_h = zeros(3,1);
                else
                   V3_h = Kp_3*norm(elbow_target_point - elbow_position_h) - Kd_3*norm(x3_dot_d_hm1);
                   x3_dot_d_h = V3_h * ((elbow_target_point - elbow_position_h) / norm(elbow_target_point - elbow_position_h));
                   x3_ddot_d_h = (x3_dot_d_h-x3_dot_d_hm1) / T;
                end
                m3 = length(x3_ddot_d_h);

                % x3_dot_d_h = -50*elbow_position_h(2);
                % x3_ddot_d_h = (x3_dot_d_h - J3_h(2,:)*q_dot_h)/T;
                % m3 = length(x3_ddot_d_h);

                % % SNS solution                
                q_dot_new = SNS_velocity_multitask(n, {m1}, {J1_h}, {x1_dot_d_h}, bounds, q_h, T, false);                                
                %q_dot_new = SNS_velocity_multitask(n, {m1, m2}, {J1_h, J2_h}, {x1_ddot_d_h, q_dot_cs}, bounds, q_h, T, false);                                                
                %q_dot_new = SNS_velocity_multitask(n, {m1, m3}, {J1_h, J3_h}, {x1_ddot_d_h, x3_ddot_d_h}, bounds, q_h, T, false);    

                q_new = q_h + q_dot_new*T;

                bounds_min_acceleration = bounds{3}(1,:);
                bounds_max_acceleration = bounds{3}(2,:);                
                for i=1:n
                    if round(q_ddot_new(i),4)>round(bounds_max_acceleration(i),4) || round(q_ddot_new(i),4)<round(bounds_min_acceleration(i),4)                        
                        saturation_counter = saturation_counter+1;
                    end
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d\n', k);
                fprintf('norm(ee_position_h - x_d) = ');disp(norm(ee_position_h-x_d))      
                fprintf('\n');
                fprintf('x_d                = ');disp(x_d');
                fprintf('ee_position_h      = ');disp(ee_position_h');                
                fprintf('elbow_target_point = ');disp(elbow_target_point')
                fprintf('elbow_position_h   =');disp(elbow_position_h');                
                fprintf('\n');
                fprintf('x1_dot_d_h         = ');disp(x1_dot_d_h');
                fprintf('x1_ddot_d_h        = ');disp(x1_ddot_d_h');   
                fprintf('x2_dot_d_h         = ');disp(x2_dot_d_h');
                fprintf('x2_ddot_d_h        = ');disp(x2_ddot_d_h');   
                fprintf('\n');
                fprintf('q_ddot_new         = ');disp(q_ddot_new')
                fprintf('q_dot_new          = ');disp(q_dot_new')
                fprintf('q_new              = ');disp(q_new')
                fprintf('saturation_counter = ');disp(saturation_counter) 
                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_dot_hm1 = q_dot_h;
                ee_position_hm1 = ee_position_h;
                J1_hm1 = J1_h;   

                J3_hm1 = J3_h;
                x3_dot_d_hm1 = x3_dot_d_h;
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
               
                J1_h = self.robot.get_J_ee(q_h);
                J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);                      

                J3_h = self.robot.get_J_eb(q_h);
                J3_dot_h = self.robot.get_J_dot_eb(q_h, q_dot_h);
                elbow_position_h = self.robot.get_elbow_position(q_h);                      
             
                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));

                % Save new configuration
                joints_positions = [joints_positions, q_h];
                joints_velocities = [joints_velocities, q_dot_h];                
                directional_errors = [directional_errors, e_d];
            end
        end

        %% run simulation acceleration level
        function [joints_positions, joints_velocities, joints_accelerations, directional_errors] = run_simulation_acceleration_level(self)
            
            n = self.robot.n;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
                            
            T = self.simulation_step;
            
            % Current configuration (time h)
            q_h = self.robot.q_0;
            q_dot_h = self.robot.q_dot_0;
            ee_position_h = self.robot.ee_position_0;            
            elbow_position_h = self.robot.elbow_position_0;

            % Jacobians of task 1, ee_position (time h)
            J1_h = self.robot.get_J_ee(q_h);
            J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);

            % Jacobians of task 2 (time h)
            J2_h = eye(n);
            J2_dot_h = zeros(7,7);

            % Jacobians of task 3, elbow_position (time h)
            J3_h = self.robot.get_J_eb(q_h);
            J3_dot_h = self.robot.get_J_dot_eb(q_h,q_dot_h);

            % Previous configuration (time h-1)
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;

            % Jacobian of task 1, ee_position (time h-1)
            J1_hm1 = self.robot.get_J_ee(q_hm1);                 
            J3_hm1 = self.robot.get_J_eb(q_hm1);   

            % Task 3:
            x3_dot_d_hm1 = zeros(3,1);      

            % Set PD controller gains of task 1
            Kp_1 = 2;
            Kd_1 = 0.05;
            
            % Set PD controller gains of task 3
            Kp_3 = 2;
            Kd_3 = 0.05;
       
            % All configurations during task
            joints_positions = [q_h];
            joints_velocities = [q_dot_h]; 
            joints_accelerations = [zeros(n,1)]; 
            directional_errors = [];
            
            k = 1;
            x_d = self.path(1:3,k);
            
            saturation_counter = 0;

            while true                
                if norm(ee_position_h - x_d) < 0.005
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                V_h = Kp_1*norm(x_d - ee_position_h) - Kd_1*norm(J1_hm1*q_dot_hm1);
                x1_dot_d_h = V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h));
                x1_ddot_d_h =  (x1_dot_d_h-J1_hm1*q_dot_hm1) / T;
                m1 = length(x1_ddot_d_h);               

                % TASK 2: self motion dumping
                q_ddot_cs = -1000*q_dot_h;
                m2 = length(q_ddot_cs);

                % TASK 3: elbow
                elbow_target_point = [elbow_position_h(1); 0; elbow_position_h(3)];                
                if(norm(elbow_target_point-elbow_position_h)==0)
                   x3_dot_d_h = zeros(3,1);
                   x3_ddot_d_h = zeros(3,1);
                else
                   V3_h = Kp_3*norm(elbow_target_point - elbow_position_h) - Kd_3*norm(x3_dot_d_hm1);
                   x3_dot_d_h = V3_h * ((elbow_target_point - elbow_position_h) / norm(elbow_target_point - elbow_position_h));
                   x3_ddot_d_h = (x3_dot_d_h-x3_dot_d_hm1) / T;
                end
                m3 = length(x3_ddot_d_h);

                % x3_dot_d_h = -50*elbow_position_h(2);
                % x3_ddot_d_h = (x3_dot_d_h - J3_h(2,:)*q_dot_h)/T;
                % m3 = length(x3_ddot_d_h);
                        
                % SNS solution
                %q_ddot_new = SNS_acceleration_multitask(n, {m1}, {J1_h}, {J1_dot_h}, {x1_ddot_d_h}, bounds, q_h, q_dot_h, T, true);                                
                q_ddot_new = SNS_acceleration_multitask(n, {m1, m2}, {J1_h, J2_h}, {J1_dot_h, J2_dot_h}, {x1_ddot_d_h, q_ddot_cs}, bounds, q_h, q_dot_h, T, false);                                                
                %q_ddot_new = SNS_acceleration_multitask(n, {m1, m3}, {J1_h, J3_h}, {J1_dot_h, J3_dot_h}, {x1_ddot_d_h, x3_ddot_d_h}, bounds, q_h, q_dot_h, T, false);                                                                                

                q_dot_new = q_dot_h + q_ddot_new*T;
                q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;

                bounds_min_acceleration = bounds{3}(1,:);
                bounds_max_acceleration = bounds{3}(2,:);                
                for i=1:n
                    if round(q_ddot_new(i),4)>round(bounds_max_acceleration(i),4) || round(q_ddot_new(i),4)<round(bounds_min_acceleration(i),4)                        
                        saturation_counter = saturation_counter+1;
                    end
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d\n', k);
                fprintf('norm(ee_position_h - x_d) = ');disp(norm(ee_position_h-x_d))      
                fprintf('\n');
                fprintf('x_d                = ');disp(x_d');
                fprintf('ee_position_h      = ');disp(ee_position_h');                
                fprintf('elbow_target_point = ');disp(elbow_target_point')
                fprintf('elbow_position_h   =');disp(elbow_position_h');                
                fprintf('\n');
                fprintf('x1_dot_d_h         = ');disp(x1_dot_d_h');
                fprintf('x1_ddot_d_h        = ');disp(x1_ddot_d_h');                  
                fprintf('x3_dot_d_h         = ');disp(x3_dot_d_h');
                fprintf('x3_ddot_d_h        = ');disp(x3_ddot_d_h');   
                fprintf('\n');
                fprintf('q_ddot_new         = ');disp(q_ddot_new')
                fprintf('q_dot_new          = ');disp(q_dot_new')
                fprintf('q_new              = ');disp(q_new')
                fprintf('saturation_counter = ');disp(saturation_counter) 
                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_dot_hm1 = q_dot_h;
                ee_position_hm1 = ee_position_h;
                J1_hm1 = J1_h;   

                J3_hm1 = J3_h;
                x3_dot_d_hm1 = x3_dot_d_h;
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
               
                J1_h = self.robot.get_J_ee(q_h);
                J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);                      

                J3_h = self.robot.get_J_eb(q_h);
                J3_dot_h = self.robot.get_J_dot_eb(q_h, q_dot_h);
                elbow_position_h = self.robot.get_elbow_position(q_h);                      
             
                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));

                % Save new configuration
                joints_positions = [joints_positions, q_h];
                joints_velocities = [joints_velocities, q_dot_h]; 
                joints_accelerations = [joints_accelerations, q_dot_h]; 
                directional_errors = [directional_errors, e_d];            
            end
        end
    end


end