classdef Simulation_Planar4R
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
        function self = Simulation_Planar4R(level, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon)

            format short
            self.robot_name = 'Planar4R';
            
            % Check parameters
            if isnan(level)
                self.level = 'velocity';
            else
                self.level = level;
            end

            
            if isnan(q_0)
                self.q_0 = [0; pi/2; -pi/2; pi/2;];                
            else
                self.q_0 = q_0;
            end

            if isnan(q_dot_0)                                                
                self.q_dot_0 = zeros(4,1);                
            else
                self.q_dot_0 = q_dot_0;    
            end

            if isnan(q_ddot_0)
                self.q_ddot_0 = zeros(4,1);                
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

            % Loading robot
            fprintf('Loading robot ... ');
            self.robot = Planar4R(self.q_0, self.q_dot_0, self.q_ddot_0);            
            fprintf('done! Robot: %s\n', self.robot_name);
                                
            % compute path
            fprintf('Creating path ... ')
            self.path = self.create_hexagonal_path(1);            
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
        
        %% Create hexagonal path
        function path = create_hexagonal_path(self, n_cycle)            
            poly_0 = [[0; 0.2; 0], [0.1732; 0.1; 0], [0.1732; -0.1; 0], [0; -0.2; 0], [-0.1732; -0.1; 0], [-0.1732; 0.1; 0]];           
            center = [0.7; 0; 0];  

            points = [];
            for n=1:n_cycle
                for i=1:6
                    points = [points, [poly_0(:,i)]];
                end
            end
            path = [points, poly_0(:,1)]+center;                       
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
            link_1_position_h = self.robot.link_1_position_0;            

            % Jacobians of task 1, ee_position (time h)
            J1_h = self.robot.get_J_ee(q_h);            

            % Jacobians of task 2, link_1_position (time h)
            J2_h = self.robot.get_J_link_1(q_h);            

            % Previous configuration (time h-1)
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;
       
            % Jacobian of task 1, ee_position (time h-1)
            J1_hm1 = self.robot.get_J_ee(q_hm1);    
            
            % Jacobian of task 2, ee_position (time h-1)
            J2_hm1 = self.robot.get_J_link_1(q_hm1); 

            % Set PD controller gains of task 1
            Kp_1 = 2;
            Kd_1 = 0.05;

            % Set PD controller gains of task 2
            Kp_2 = 2;
            Kd_2 = 0.05;
                 
            % All configurations during task
            joints_positions = [q_h];  
            joints_velocities = [q_dot_h];            
            directional_errors = [];
            
            k = 1;
            x_d = self.path(1:3,k);
            link_1_target_position = [0.25;0;0];
            
            saturation_counter = 0;

            while true                
                if norm(ee_position_h - x_d) < self.epsilon
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                V1_h = Kp_1*norm(x_d - ee_position_h) - Kd_1*norm(J1_hm1*q_dot_hm1);
                x1_dot_d_h = V1_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h));
                x1_ddot_d_h =  (x1_dot_d_h-J1_hm1*q_dot_hm1) / T;
                m1 = length(x1_ddot_d_h);               
 
                % TASK 2                
                V2_h = Kp_2*norm(link_1_target_position - link_1_position_h) - Kd_2*norm(J2_hm1*q_dot_hm1);
                x2_dot_d_h = V2_h * ((link_1_target_position - link_1_position_h) / norm(link_1_target_position - link_1_position_h));
                x2_ddot_d_h =  (x2_dot_d_h-J2_hm1*q_dot_hm1) / T;
                m2 = length(x2_ddot_d_h); 

                % SNS solution single task (first task)
                q_dot_new = SNS_velocity_multitask(n, {m1}, {J1_h}, {x1_dot_d_h}, bounds, q_h, T, false);      
                
                % SNS solution multiple task (first + second task)
                %q_dot_new = SNS_velocity_multitask(n, {m1, m2}, {J1_h, J2_h}, {x1_dot_d_h, x2_dot_d_h}, bounds, q_h, T, false); 

                q_new = q_h + q_dot_new*T;

                bounds_min_velocity = bounds{3}(1,:);
                bounds_max_velocity = bounds{3}(2,:);               
                for i=1:n
                    if round(q_dot_new(i),4)>round(bounds_max_velocity(i),4) || round(q_dot_new(i),4)<round(bounds_min_velocity(i),4)                        
                        saturation_counter = saturation_counter+1;
                    end
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d\n', k);
                fprintf('norm(ee_position_h - x_d) = ');disp(norm(ee_position_h-x_d))      
                fprintf('\n');
                fprintf('x_d               = ');disp(x_d');
                fprintf('ee_position_h     = ');disp(ee_position_h');
                fprintf('link_1_position_h = ');disp(link_1_position_h');              
                fprintf('\n');
                fprintf('x1_dot_d_h        = ');disp(x1_dot_d_h');
                fprintf('x1_ddot_d_h       = ');disp(x1_ddot_d_h');   
                fprintf('x2_dot_d_h        = ');disp(x2_dot_d_h');
                fprintf('x2_ddot_d_h       = ');disp(x2_ddot_d_h');   
                fprintf('\n');
                fprintf('q_ddot_new        = ');disp(q_ddot_new')
                fprintf('q_dot_new         = ');disp(q_dot_new')
                fprintf('q_new             = ');disp(q_new')
                fprintf('saturation_counter= ');disp(saturation_counter) 
                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_dot_hm1 = q_dot_h;                
                
                J1_hm1 = J1_h;                   
                J2_hm1 = J2_h;   
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
               
                J1_h = self.robot.get_J_ee(q_h);                
                ee_position_h = self.robot.get_ee_position(q_h);                      
                
                J2_h = self.robot.get_J_link_1(q_h);                
                link_1_position_h = self.robot.get_link_1_position(q_h);  
                                   
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
            link_1_position_h = self.robot.link_1_position_0;            

            % Jacobians of task 1, ee_position (time h)
            J1_h = self.robot.get_J_ee(q_h);
            J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);

            % Jacobians of task 2, link_1_position (time h)
            J2_h = self.robot.get_J_link_1(q_h);
            J2_dot_h = self.robot.get_J_dot_link_1(q_h, q_dot_h);

            % Previous configuration (time h-1)
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;
    
            % Jacobian of task 1, ee_position (time h-1)
            J1_hm1 = self.robot.get_J_ee(q_hm1);    
            
            % Jacobian of task 2, ee_position (time h-1)
            J2_hm1 = self.robot.get_J_link_1(q_hm1); 

            % Set PD controller gains of task 1
            Kp_1 = 2;
            Kd_1 = 0.05;

            % Set PD controller gains of task 2
            Kp_2 = 2;
            Kd_2 = 0.05;
                 
            % All configurations during task
            joints_positions = [q_h];  
            joints_velocities = [q_dot_h];
            joints_accelerations = [zeros(n,1)];
            directional_errors = [];
            
            k = 1;
            x_d = self.path(1:3,k);
            link_1_target_position = [0.25;0;0];
            
            saturation_counter = 0;

            while true                
                if norm(ee_position_h - x_d) < self.epsilon
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                V1_h = Kp_1*norm(x_d - ee_position_h) - Kd_1*norm(J1_hm1*q_dot_hm1);
                x1_dot_d_h = V1_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h));
                x1_ddot_d_h =  (x1_dot_d_h-J1_hm1*q_dot_hm1) / T;
                m1 = length(x1_ddot_d_h);               
 
                % TASK 2                
                V2_h = Kp_2*norm(link_1_target_position - link_1_position_h) - Kd_2*norm(J2_hm1*q_dot_hm1);
                x2_dot_d_h = V2_h * ((link_1_target_position - link_1_position_h) / norm(link_1_target_position - link_1_position_h));
                x2_ddot_d_h =  (x2_dot_d_h-J2_hm1*q_dot_hm1) / T;
                m2 = length(x2_ddot_d_h);            


                % SNS solution single task (first task)
                %q_ddot_new = SNS_acceleration_multitask(n, {m1}, {J1_h}, {J1_dot_h}, {x1_ddot_d_h}, bounds, q_h, q_dot_h, T, false);                                
                
                % SNS solution multiple task (first + second task)
                q_ddot_new = SNS_acceleration_multitask(n, {m1, m2}, {J1_h, J2_h}, {J1_dot_h, J2_dot_h}, {x1_ddot_d_h, x2_ddot_d_h}, bounds, q_h, q_dot_h, T, false);                                                
                
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
                fprintf('x_d               = ');disp(x_d');
                fprintf('ee_position_h     = ');disp(ee_position_h');
                fprintf('link_1_position_h = ');disp(link_1_position_h');              
                fprintf('\n');
                fprintf('x1_dot_d_h        = ');disp(x1_dot_d_h');
                fprintf('x1_ddot_d_h       = ');disp(x1_ddot_d_h');   
                fprintf('x2_dot_d_h        = ');disp(x2_dot_d_h');
                fprintf('x2_ddot_d_h       = ');disp(x2_ddot_d_h');   
                fprintf('\n');
                fprintf('q_ddot_new        = ');disp(q_ddot_new')
                fprintf('q_dot_new         = ');disp(q_dot_new')
                fprintf('q_new             = ');disp(q_new')
                fprintf('saturation_counter= ');disp(saturation_counter) 
                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_dot_hm1 = q_dot_h;
                
                J1_hm1 = J1_h;   
                J2_hm1 = J2_h;   
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
                q_ddot_h = q_ddot_new;
               
                J1_h = self.robot.get_J_ee(q_h);
                J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);                      
                
                J2_h = self.robot.get_J_link_1(q_h);
                J2_dot_h = self.robot.get_J_dot_link_1(q_h, q_dot_h);
                link_1_position_h = self.robot.get_link_1_position(q_h);  
                                   
                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));

                % Save new configuration
                joints_positions = [joints_positions, q_h];
                joints_velocities = [joints_velocities, q_dot_h];
                joints_accelerations = [joints_accelerations, q_ddot_h];
                directional_errors = [directional_errors, e_d];
            
            end
        end
    end
end