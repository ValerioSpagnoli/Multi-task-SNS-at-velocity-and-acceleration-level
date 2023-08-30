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
        function self = Simulation(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step)

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
                self.q_0 = [0; pi/4; pi/4; pi/4; 0; 0; 0];
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
                     
            fprintf('Loading robot ... ');
            if strcmp(self.robot_name, 'KUKA_LBR_IIWA_7_R800')
                self.robot = KUKA_LBR_IIWA7(self.q_0, self.q_dot_0, self.q_ddot_0);
            elseif strcmp(self.robot_name, 'KUKA_LBR_IV')
                self.robot = KUKA_LBR_IV(self.q_0, self.q_dot_0, self.q_ddot_0);
            end
            fprintf('done! \n');
            
            fprintf('Creating path ... ')
            self.path = self.create_path(1);
            fprintf('done! \n');
            
            fprintf('Start simulation ... \n\n')
            pause(2);
            if strcmp(level, 'velocity')
                [self.joint_positions, self.directional_error] = self.run_simulation_velocity_level();
            elseif strcmp(level, 'acceleration')
                [self.joint_positions, self.directional_error] = self.run_simulation_acceleration_level();
            end
            fprintf('Simulation ended. \n')
        end

        %% create path
        function path = create_path(self, n_cycle)
            poly_0 = [[0; 0; 0.2], [0; 0.1732; 0.1], [0; 0.1732; -0.1], [0; 0; -0.2], [0; -0.1732; -0.1], [0; -0.1732; 0.1]];
            center = [0.1; 0.35; 0.6235];
            
            points = [];
            for n=1:n_cycle
                for i=1:6
                    points = [points, [poly_0(:,i)]];
                end
            end
            path = [points, poly_0(:,1)]+center;
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
            q_ddot_h = self.q_ddot_0;
            ee_position_h = self.robot.ee_position_0;
            J1_h = self.robot.get_J(self.q_0);
            
            % configuration at time h-1
            q_hm1 = self.q_0;
            q_dot_hm1 = self.q_dot_0;
            q_ddot_hm1 = self.q_ddot_0;
            ee_position_hm1 = self.robot.ee_position_0;
            J1_hm1 = self.robot.get_J(self.q_0);
            
            % all configurations during task
            qs = [q_h];
            eds = [];
        
            k = 1;
            while k<=size(self.path,2)
                
                x_d = round(self.path(1:3,k),4);
                
                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n', k, norm(ee_position_h-x_d));
            
                if norm(ee_position_h - x_d) < 0.01
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

                bounds_min_velocity = bounds{2}(1,:);
                bounds_max_velocity = bounds{2}(2,:);
                for i=1:ndof
                    if round(q_dot_new(i),4)>round(bounds_max_velocity(i),4)
                        fprintf('Joint %d out of max velocity bounds: q_dot(%d) = %f, bounds_max_velocity(%d) = %f\n',i, i, q_dot_new(i), i, bounds_max_velocity(i));
                    elseif round(q_dot_new(i),4)<round(bounds_min_velocity(i),4)
                        fprintf('Joint %d out of min velocity bounds: q_dot(%d) = %f, bounds_min_velocity(%d) = %f\n',i, i, q_dot_new(i), i, bounds_min_velocity(i));
                    end
                end
        
                % set previous configuration to current configuration ...
                q_hm1 = q_h;
                q_dot_hm1 = q_dot_h;
                J1_hm1 = J1_h;
        
                % update current cunfiguration to new configuration
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

            bounds_min_position = bounds{1}(1,:);
            bounds_max_position = bounds{1}(2,:);
            
            bounds_min_velocity = bounds{2}(1,:);
            bounds_max_velocity = bounds{2}(2,:);
                    
            bounds_min_acceleration = bounds{3}(1,:);
            bounds_max_acceleration = bounds{3}(2,:);

                            
            T = self.simulation_step;
            kp = 10;
            kd = 0.1;
            
            % current configuration (time h)
            q_h = self.robot.q_0;
            q_dot_h = self.robot.q_dot_0;
            q_ddot_h = self.robot.q_ddot_0;
            ee_position_h = self.robot.ee_position_0;

            J1_h = self.robot.get_J(q_h);
            J1_dot_h = self.robot.get_J_dot(q_h, q_dot_h);

            J2_h = eye(ndof);
            J2_dot_h = zeros(7,7);

            % configuration at time h-1
            q_hm1 = self.robot.q_0;
            q_dot_hm1 = self.robot.q_ddot_0;
            q_ddot_hm1 = self.robot.q_ddot_0;
            ee_position_hm1 = self.robot.ee_position_0;

            J1_hm1 = self.robot.get_J(q_hm1);
            J1_dot_hm1 = self.robot.get_J_dot(q_hm1, q_dot_hm1);

            J2_hm1 = eye(ndof);
            J2_dot_hm1 = zeros(7,7);
            
            % all configurations during task
            qs = [q_h];
            eds = [];
        
            k = 1;
            while k<=size(self.path,2)
                
                x_d = round(self.path(1:3,k),4);
                    
                if norm(ee_position_h - x_d) < 0.01
                    k = k+1;
                    x_d = self.path(k);            
                end
        
                % TASK 1: path following
                V_h = round(kp*norm(x_d - ee_position_h) - kd*norm(J1_hm1*q_dot_hm1),4);
                x_dot_d_h = round(V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h)),4);
                x_ddot_d_h = round( (x_dot_d_h-J1_hm1*q_dot_hm1) / T ,4);
                
                % A_h = round(kp*norm(x_dot_d_h - J_h*q_dot_h) - kd*norm(J_hm1*q_ddot_hm1+J_dot_hm1*q_dot_hm1),4);
                % x_ddot_d_h = round(A_h * ((x_dot_d_h - J_h*q_dot_h) / norm(x_dot_d_h - J_h*q_dot_h)),4);

                m1 = length(x_ddot_d_h);


                % TASK 2: self motion dumping
                q_dot_cs = -1000*q_dot_h;

                m2 = length(q_dot_cs);


                % Shaping joint bounds
                [bounds_Q_ddot_min, bounds_Q_ddot_max] = self.shaping_acceleration_bounds(ndof, bounds, q_h, q_dot_h, T);
                bounds = {[bounds_min_position; bounds_max_position], [bounds_min_velocity; bounds_max_velocity], [bounds_Q_ddot_min; bounds_Q_ddot_max]};


                % SNS solution
                q_ddot_new = SNS_acceleration_multitask_full(ndof, {m1}, {J1_h}, {J1_dot_h}, {x_ddot_d_h}, bounds, q_h, q_dot_h, T, true);
                q_dot_new = q_dot_h + q_ddot_new*T;
                q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;


                % Pseudoinversion solution
                % q_ddot_new = pseudoinverse_solution(m, ndof, {J_h}, {J_dot_h}, {x_ddot_d_h}, bounds, q_dot_h, true);
                % q_dot_new = q_dot_h + q_ddot_new*T;
                % q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;

                
                % hard bounds
                limit_exceeded_acc = false;
                limit_exceeded_vel = false;
                limit_exceeded_pos = false;

                for i=1:ndof
                    if round(q_ddot_new(i),4)>round(bounds_max_acceleration(i),4)
                        limit_exceeded_acc=true;
                    elseif round(q_ddot_new(i),4)<round(bounds_min_acceleration(i),4)
                        limit_exceeded_acc=true;
                    end
                    
                    if round(q_dot_new(i),4)>round(bounds_max_velocity(i),4)
                        %q_dot_new(i) = round(bounds_max_velocity(i),4);
                        limit_exceeded_vel=true;
                    elseif round(q_dot_new(i),4)<round(bounds_min_velocity(i),4)
                        %q_dot_new(i) = round(bounds_min_velocity(i),4);
                        limit_exceeded_vel=true;
                    end
                    
                    if round(q_new(i),4)>round(bounds_max_position(i),4)
                        %q_new(i) = round(bounds_max_position(i),4);
                        limit_exceeded_pos=true;
                    elseif round(q_new(i),4)<round(bounds_min_position(i),4)
                        %q_new(i) = round(bounds_min_position(i),4);
                        limit_exceeded_pos=true;
                    end
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n\n', k, norm(ee_position_h-x_d));

                fprintf('ee_position_h = ');disp(ee_position_h');
                fprintf('q_dot_hm1 = ');disp(q_dot_hm1');
                fprintf('J_hm1 = \n');disp(J1_hm1);
                fprintf('J_hm1*q_dot_hm1 = ');disp(round(J1_hm1*q_dot_hm1,4)');

                fprintf('V_h = ');disp(V_h);
                fprintf('x_d = ');disp(x_d');
                fprintf('x_dot_d_h = ');disp(x_dot_d_h');
                fprintf('x_ddot_d_h = ');disp(x_ddot_d_h');
                
                fprintf('q_ddot_SNS = ');disp(q_ddot_new')
                fprintf('q_dot_SNS = ');disp(q_dot_new')
                fprintf('q_SNS = ');disp(q_new')

                fprintf('Limit_exceeded_acc = ');disp(limit_exceeded_acc);
                fprintf('Limit_exceeded_vel = ');disp(limit_exceeded_vel);
                fprintf('Limit_exceeded_pos = ');disp(limit_exceeded_pos);

                fprintf('==============================================================================\n')

                if isnan(q_ddot_new)
                    break;
                end

                % set previous configuration to current configuration ...
                q_hm1 = q_h;
                q_dot_hm1 = q_dot_h;
                q_ddot_hm1 = q_ddot_h;
                J1_hm1 = J1_h;
                J1_dot_hm1 = J1_dot_h;
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
                q_ddot_h = q_ddot_new;

                J1_h = self.robot.get_J(q_h);
                J1_dot_h = self.robot.get_J_dot(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);

                % Directional error
                e_d = acos(((x_d-ee_position_h)/norm(x_d-ee_position_h)) .* ((J1_h*q_dot_h)/norm(J1_h*q_dot_h)));
        
                % save new configuration
                qs = [qs, q_h];
                eds = [eds, e_d];
            end
        end

        %% Shaping joint bounds
        function [bounds_Q_ddot_min, bounds_Q_ddot_max] = shaping_acceleration_bounds(self, n, bounds, q, q_dot, T)

            bounds_min_position = bounds{1}(1,:);
            bounds_max_position = bounds{1}(2,:);
        
            bounds_min_velocity = bounds{2}(1,:);
            bounds_max_velocity = bounds{2}(2,:);
        
            bounds_min_acceleration = bounds{3}(1,:);
            bounds_max_acceleration = bounds{3}(2,:);
        
            bounds_Q_ddot_min = zeros(1,7);
            bounds_Q_ddot_max = zeros(1,7);
        
            for i=1:n
        
                Q_ddot_min_1 = 2*(bounds_min_position(i)-q(i)-q_dot(i)*T)/T^2;
                Q_ddot_min_2 = -(bounds_max_velocity(i)+q_dot(i))/T;
                Q_ddot_min_3 = bounds_min_acceleration(i);

                if Q_ddot_min_1 > 0
                    Q_ddot_min_1 = -Inf;
                end
                if Q_ddot_min_2 > 0
                    Q_ddot_min_2 = -Inf;
                end
                
                bounds_Q_ddot_min(i) = round( max([Q_ddot_min_1, Q_ddot_min_2, Q_ddot_min_3]) ,4);
        
                
                Q_ddot_max_1 = 2*(bounds_max_position(i)-q(i)-q_dot(i)*T)/T^2;
                Q_ddot_max_2 = (bounds_max_velocity(i)-q_dot(i))/T;
                Q_ddot_max_3 = bounds_max_acceleration(i);

                if Q_ddot_max_1 < 0
                    Q_ddot_max_1 = +Inf;
                end
                if Q_ddot_max_2 > 0
                    Q_ddot_max_2 = +Inf;
                end
                
                bounds_Q_ddot_max(i) = round( min([Q_ddot_max_1, Q_ddot_max_2, Q_ddot_max_3]) ,4);   

            end    
        end
    end


end