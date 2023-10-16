classdef Simulation
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
        round_up
        round_point
        joints_positions
        directional_errors
        elbow_positions
        elbow_velocities
    end

    methods
        %% Constructur
        % level = {'acceleration', 'velocity'}
        % robot_name = {'KUKA_LBR_IIWA_7_R800', 'KUKA_LBR_IV'}
        % q_0: [7x1] vector [rad]
        % q_dot_0: [7x1] vector [rad/s]
        % q_ddot_0: [7x1] vector [rad/s^2]
        % simulation_step: time step of simulation
        function self = Simulation(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon, round_up, round_point)

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
                self.q_dot_0 = zeros(7,1);
            else
                self.q_dot_0 = q_dot_0;    
            end
            if isnan(q_ddot_0)
                self.q_ddot_0 = zeros(7,1);
            else
                self.q_ddot_0 = q_ddot_0;
            end
            if isnan(simulation_step)
                self.simulation_step = 0.001;
            else
                self.simulation_step = simulation_step;
            end
            if isnan(epsilon)
                self.epsilon = 0.01;
            else
                self.epsilon = epsilon;
            end
            if isnan(round_up)
                self.round_up = false;
            else
                self.round_up = round_up;
            end
            if isnan(round_point)
                self.round_point = 10;
            else
                self.round_point = round_point;
            end
                     
            % load robot
            fprintf('Loading robot ... ');
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                self.robot = KUKA_LBR_IIWA7(self.q_0, self.q_dot_0, self.q_ddot_0);
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                self.robot = KUKA_LBR_IV(self.q_0, self.q_dot_0, self.q_ddot_0);
            end
            fprintf('done! Robot: %s\n', robot_name);
            
            % compute path
            fprintf('Creating path ... ')
            self.path = self.create_hexagonal_path(1);
            %self.path = self.create_circular_path(1);
            fprintf('done! \n');

           
            % compute simulation
            fprintf('Start simulation ... \n\n')
            pause(1);
            if strcmp(level, 'velocity')
                [self.joints_positions, self.directional_errors] = self.run_simulation_velocity_level();
            elseif strcmp(level, 'acceleration')
                [self.joints_positions, self.directional_errors, self.elbow_positions, self.elbow_velocities] = self.run_simulation_acceleration_level();
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
                    points = [points, [poly_0(:,i)]];
                end
            end
            path = [points, poly_0(:,1)]+center;
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
        function [qs, eds] = run_simulation_velocity_level(self)

            ndof = self.robot.ndof;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
                            
            T = self.simulation_step;
            Kp = 10;
            Kd = 0.1;
            
            % current configuration (time h)
            q_h = self.q_0;
            q_dot_h = self.q_dot_0;

            ee_position_h = self.robot.ee_position_0;
            J1_h = self.robot.get_J_ee(self.q_0);

            elbow_position_h = self.robot.elbow_position_0;
            J2_h = self.robot.get_J_eb(self.q_0);
            
            % configuration at time h-1
            q_dot_hm1 = self.q_dot_0;
            ee_position_hm1 = zeros(3,1);
            J1_hm1 = self.robot.get_J_ee(self.q_0);
            J2_hm1 = self.robot.get_J_eb(self.q_0);
            
            % all configurations during task
            qs = [q_h];
            eds = [];
        
            k = 1;
            x_d = round(self.path(1:3,k),4);

            while true

                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n', k, norm(ee_position_h-x_d));
            
                if norm(ee_position_h - x_d) < 0.005
                    k = k+1;
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = round(self.path(1:3,k),4);
                end
        
                % TASK 1: path_following
                V = round(Kp*norm(x_d-ee_position_h) - Kd*norm(J1_hm1*q_dot_hm1),4);
                x1_dot_d_h = round(V * ((x_d-ee_position_h)/norm(x_d-ee_position_h)),4);
                m1 = length(x_d);

                % TASK 2: elbow position close to x-z plane
                V_h = round( Kp*norm(elbow_position_h(2,:)) - Kd*norm(J2_hm1(2,:)*q_dot_hm1) ,4);
                x2_dot_d_h = round( V_h * (elbow_position_h(2,:) / norm(elbow_position_h(2,:))) ,4);
                m2 = length(x2_dot_d_h);
                % x2_dot_d_h = round( -50 * elbow_position_h(2,:), 4);
                % m2 = length(x2_dot_d_h);

                % % SNS solution
                %q_dot_new = SNS_velocity_multitask(ndof, {m1, m2}, {J1_h, J2_h(2,:)}, {x1_dot_d_h, x2_dot_d_h}, bounds, q_h, T, false);                
                q_dot_new = SNS_velocity_multitask(ndof, {m1}, {J1_h}, {x1_dot_d_h}, bounds, q_h, T, false);                
                q_new = q_h + q_dot_new*T;

                % Set previous configuration to current configuration
                q_dot_hm1 = q_dot_h;
                J1_hm1 = J1_h;
                J2_hm1 = J2_h;
        
                % Update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
                ee_position_h = self.robot.get_ee_position(q_h);
                J1_h = self.robot.get_J_ee(q_h);
                elbow_position_h = self.robot.get_elbow_position(q_h);
                J2_h = self.robot.get_J_eb(q_h);

                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));
        
                % save new configuration
                qs = [qs, q_h];
                eds = [eds, e_d];

            end
        end

        %% run simulation acceleration level
        function [joints_positions, directional_errors, elbow_positions, elbow_velocities] = run_simulation_acceleration_level(self)

            ndof = self.robot.ndof;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
                            
            T = self.simulation_step;
            Kp = 3;
            Kd = 0.05;
            
            % Current configuration (time h)
            q_h = self.robot.q_0;
            q_dot_h = self.robot.q_dot_0;
            ee_position_h = self.robot.ee_position_0;
            elbow_position_h = self.robot.elbow_position_0;

            % Jacobians of task 1, ee_position (time h)
            J1_h = self.robot.get_J_ee(q_h);
            J1_dot_h = self.robot.get_J_dot_ee(q_h, q_dot_h);

            % Jacobians of task 2 (time h)
            J2_h = eye(ndof);
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
       
            % All configurations during task
            joints_positions = [q_h];
            directional_errors = [];
            elbow_positions = [elbow_position_h];
            elbow_velocities = [q_dot_h(4)];
            
            k = 1;
            x_d = self.path(1:3,k);
            
            counter = 0;

            while true

                if norm(ee_position_h - x_d) < 0.005
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                if self.round_up
                V_h = round( Kp*norm(x_d - ee_position_h) - Kd*norm(J1_hm1*q_dot_hm1) ,self.round_point);
                x1_dot_d_h = round( V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h)) ,self.round_point);
                x1_ddot_d_h =  round( (x1_dot_d_h-J1_hm1*q_dot_hm1) / T ,self.round_point);
                m1 = length(x1_ddot_d_h);
                else
                V_h = Kp*norm(x_d - ee_position_h) - Kd*norm(J1_hm1*q_dot_hm1);
                x1_dot_d_h = V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h));
                x1_ddot_d_h =  (x1_dot_d_h-J1_hm1*q_dot_hm1) / T;
                m1 = length(x1_ddot_d_h);
                end

                % TASK 2: self motion dumping
                q_ddot_cs = -1000*q_dot_h;
                m2 = length(q_ddot_cs);

                % TASK 3: elbow
                elbow_target_point = [elbow_position_h(1); 0; elbow_position_h(3)];                
                if(norm(elbow_target_point-elbow_position_h)==0)
                    x3_dot_d_h = zeros (3,1);
                    x3_ddot_d_h = zeros (3,1);
                else
                    if self.round_up
                    V3_h = round( Kp*norm(elbow_target_point - elbow_position_h) - Kd*norm(J3_hm1*q_dot_hm1) ,self.round_point);
                    x3_dot_d_h = round( V3_h * ((elbow_target_point - elbow_position_h) / norm(elbow_target_point - elbow_position_h)) ,self.round_point);
                    x3_ddot_d_h = round( (x3_dot_d_h-J3_hm1*q_dot_hm1) / T ,self.round_point);
                    else
                    V3_h = Kp*norm(elbow_target_point - elbow_position_h) - Kd*norm(x3_dot_d_hm1);
                    x3_dot_d_h = V3_h * ((elbow_target_point - elbow_position_h) / norm(elbow_target_point - elbow_position_h));
                    x3_ddot_d_h = (x3_dot_d_h-x3_dot_d_hm1) / T;
                    end
                end
                m3 = length(x3_ddot_d_h);
        
                % SNS solution
                %q_ddot_new = SNS_acceleration_multitask(ndof, {m1}, {J1_h}, {J1_dot_h}, {x1_ddot_d_h}, bounds, q_h, q_dot_h, T, self.round_up, self.round_point, true);                                
                q_ddot_new = SNS_acceleration_multitask(ndof, {m1, m2}, {J1_h, J2_h}, {J1_dot_h, J2_dot_h}, {x1_ddot_d_h, q_ddot_cs}, bounds, q_h, q_dot_h, T, self.round_up, self.round_point, false);                                                
                %q_ddot_new = SNS_acceleration_multitask(ndof, {m1, m3}, {J1_h, J3_h}, {J1_dot_h, J3_dot_h}, {x1_ddot_d_h, x3_ddot_d_h}, bounds, q_h, q_dot_h, T, self.round_up, self.round_point, false);                                                                

                q_dot_new = q_dot_h + q_ddot_new*T;
                q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;

                bounds_min_acceleration = bounds{3}(1,:);
                bounds_max_acceleration = bounds{3}(2,:);
                
                for i=1:ndof
                    if round(q_ddot_new(i),4)>round(bounds_max_acceleration(i),4) || round(q_ddot_new(i),4)<round(bounds_min_acceleration(i),4)                        
                        counter = counter+1;
                    end
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d\n', k);
                fprintf('norm(ee_position_h - x_d) = ');disp(norm(ee_position_h-x_d))      
                fprintf('\n');
                fprintf('x_d           = ');disp(x_d');
                fprintf('ee_position_h = ');disp(ee_position_h');           
                fprintf('\n');
                fprintf('elbow_target_point = ');disp(elbow_target_point');
                fprintf('elbow_position_h   =');disp(elbow_position_h');
                fprintf('\n');
                fprintf('x1_dot_d_h    = ');disp(x1_dot_d_h');
                fprintf('x1_ddot_d_h   = ');disp(x1_ddot_d_h');   
                fprintf('x3_dot_d_h    = ');disp(x3_dot_d_h');
                fprintf('x3_ddot_d_h   = ');disp(x3_ddot_d_h');   
                fprintf('\n');
                fprintf('q_ddot_new    = ');disp(q_ddot_new')
                fprintf('q_dot_new     = ');disp(q_dot_new')
                fprintf('q_new         = ');disp(q_new')

                fprintf('counter = ');disp(counter)

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
                directional_errors = [directional_errors, e_d];
                elbow_positions = [elbow_positions, elbow_position_h];
                elbow_velocities = [elbow_velocities, q_dot_h(4)];
            
            end
        end
    end


end