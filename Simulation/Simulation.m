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
        joint_positions
        directional_error
    end

    methods
        %% Constructur
        % level = {'acceleration', 'velocity'}
        % robot_name = {'KUKA_LBR_IIWA_7_R800', 'KUKA_LBR_IV'}
        % q_0: [7x1] vector [rad]
        % q_dot_0: [7x1] vector [rad/s]
        % q_ddot_0: [7x1] vector [rad/s^2]
        % simulation_step: time step of simulation
        function self = Simulation(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon)

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
                     
            % load robot
            fprintf('Loading robot ... ');
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                self.robot = KUKA_LBR_IIWA7(self.q_0, self.q_dot_0, self.q_ddot_0);
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                self.robot = KUKA_LBR_IV(self.q_0, self.q_dot_0, self.q_ddot_0);
            end
            fprintf('done! \n');
            
            % compute path
            fprintf('Creating path ... ')
            self.path = self.create_hexagonal_path(1);
            %self.path = self.create_circular_path(1);
            fprintf('done! \n');
            
            % compute simulation
            fprintf('Start simulation ... \n\n')
            pause(1);
            if strcmp(level, 'velocity')
                [self.joint_positions, self.directional_error] = self.run_simulation_velocity_level();
            elseif strcmp(level, 'acceleration')
                [self.joint_positions, self.directional_error] = self.run_simulation_acceleration_level();
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
            kp = 10;
            kd = 0.1;
            
            % current configuration (time h)
            q_h = self.q_0;
            q_dot_h = self.q_dot_0;
            ee_position_h = self.robot.ee_position_0;
            J1_h = self.robot.get_J(self.q_0);
            
            % configuration at time h-1
            q_dot_hm1 = self.q_dot_0;
            J1_hm1 = self.robot.get_J(self.q_0);
            
            % all configurations during task
            qs = [q_h];
            eds = [];
        
            k = 1;
            x_d = round(self.path(1:3,k),4);

            while k<=size(self.path,2)

                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n', k, norm(ee_position_h-x_d));
            
                if norm(ee_position_h - x_d) < 0.005
                    k = k+1;
                    x_d = self.path(k);            
                end
        
                % TASK 1: path_following
                V = round(kp*norm(x_d-ee_position_h) - kd*norm(J1_hm1*q_dot_hm1),4);
                x_dot_d_h = round(V * ((x_d-ee_position_h)/norm(x_d-ee_position_h)),4);
                m1 = length(x_d);

                % % SNS solution
                q_dot_new = SNS_velocity_multitask(ndof, {m1}, {J1_h}, {x_dot_d_h}, bounds, q_h, T, false);
                q_new = q_h + q_dot_new*T;

                % Set previous configuration to current configuration
                q_dot_hm1 = q_dot_h;
                J1_hm1 = J1_h;
        
                % Update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
                J1_h = self.robot.get_J(q_h);
                ee_position_h = self.robot.get_ee_position(q_h);

                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));
        
                % save new configuration
                qs = [qs, q_h];
                eds = [eds, e_d];

            end
        end

        %% run simulation acceleration level
        function [qs, eds] = run_simulation_acceleration_level(self)

            ndof = self.robot.ndof;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
                            
            T = self.simulation_step;
            Kp = 5;
            Kd = 0.1;
            
            % Current configuration (time h)
            q_h = self.robot.q_0;
            q_dot_h = self.robot.q_dot_0;
            q_ddot_h = self.robot.q_ddot_0;
            ee_position_h = self.robot.ee_position_0;

            % Jacobians of task 1 (time h)
            J1_h = self.robot.get_J(q_h);
            J1_dot_h = self.robot.get_J_dot(q_h, q_dot_h);

            % Jacobians of task 2 (time h)
            J2_h = eye(ndof);
            J2_dot_h = zeros(7,7);

            % Previous configuration (time h-1)
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;

            % Jacobian of task 1 (time h-1)
            J1_hm1 = self.robot.get_J(q_hm1);            

            % All configurations during task
            qs = [q_h];
            eds = [];
            
            k = 1;
            x_d = self.path(1:3,k);
            
            while true

                if norm(ee_position_h - x_d) < 0.01
                    k = k+1;       
                    if k>size(self.path,2)
                        break;
                    end
                    x_d = self.path(1:3,k);
                end
                
                % TASK 1: path following, paper formulation
                V_h = round(Kp*norm(x_d - ee_position_h) - Kd*norm(J1_hm1*q_dot_hm1),4);
                x_dot_d_h = round(V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h)),4);
                x_ddot_d_h =  round( (x_dot_d_h-J1_hm1*q_dot_hm1) / T ,4);
                m1 = length(x_ddot_d_h);

                % TASK 2: self motion dumping
                q_ddot_cs = -1000*q_dot_h;
                m2 = length(q_ddot_cs);
                
                % SNS solution
                q_ddot_new = SNS_acceleration_multitask(ndof, {m1, m2}, {J1_h, J2_h}, {J1_dot_h, J2_dot_h}, {x_ddot_d_h, q_ddot_cs}, bounds, q_h, q_dot_h, T, false);
                q_dot_new = q_dot_h + q_ddot_new*T;
                q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;

                fprintf('==============================================================================\n')

                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n\n', k, norm(ee_position_h-x_d));

                fprintf('x_d           = ');disp(x_d');
                fprintf('x_dot_d_h     = ');disp(x_dot_d_h');
                fprintf('x_ddot_d_h    = ');disp(x_ddot_d_h');
                
                fprintf('q_ddot_new    = ');disp(q_ddot_new')
                fprintf('q_dot_new     = ');disp(q_dot_new')
                fprintf('q_new         = ');disp(q_new')

                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_dot_hm1 = q_dot_h;
                J1_hm1 = J1_h;
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
               
                J1_h = self.robot.get_J(q_h);
                J1_dot_h = self.robot.get_J_dot(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);

                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));

                % Save new configuration
                qs = [qs, q_h];
                eds = [eds, e_d];
            end
        end
    end


end