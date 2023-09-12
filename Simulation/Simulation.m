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
        x_dot_d
        x_ddot_d
        versors
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
            self.path = self.create_path(1);
            fprintf('done! \n');
            
            % compute simulation
            fprintf('Start simulation ... \n\n')
            pause(1);
            if strcmp(level, 'velocity')
                [self.joint_positions, self.directional_error] = self.run_simulation_velocity_level();
            elseif strcmp(level, 'acceleration')
                [self.joint_positions, self.directional_error, self.x_dot_d, self.x_ddot_d, self.versors] = self.run_simulation_acceleration_level();
            end
            fprintf('Simulation ended. \n')
        end

        %% create path
        function path = create_path(self, n_cycle)
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
        function [qs, eds, x_dot_d, x_ddot_d, versors] = run_simulation_acceleration_level(self)

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
            
            t=0;
            e_h = 0;
            e_tot = 0;

            % all configurations during task
            qs = [q_h];
            eds = [];

            x_dot_d = [];
            x_ddot_d = [];
            versors = {};
            
            k = 1;
            x_d = self.path(1:3,k);
            pi = self.robot.ee_position_0;
            pf = x_d;
            versors{end+1} = [pi, (pf-pi)/norm(pf-pi)];

            saturation_counter = 0;
            
            while true
        
                % TASK 1: path following, paper formulation
                V_h = round(kp*norm(x_d - ee_position_h) - kd*norm(J1_hm1*q_dot_hm1),4);
                x_dot_d_h = round(V_h * ((x_d - ee_position_h) / norm(x_d - ee_position_h)),4);
                x_ddot_d_h =  round( (x_dot_d_h-J1_hm1*q_dot_hm1) / T ,4);

                m1 = length(x_ddot_d_h);
                 
                % TASK 2: self motion dumping
                % q_dot_cs = -1000*q_dot_h;
                % m2 = length(q_dot_cs);

                v_max = 0.5; 
                a_max = 2;
                %v_max = 0.4*bounds_max_velocity(2)+0.4*bounds_min_velocity(4)+0.126*bounds_max_velocity(6);
                %a_max = 0.4*bounds_max_acceleration(2)+0.4*bounds_min_acceleration(4)+0.126*bounds_max_acceleration(6);

                L = norm(pf-pi);
                v = (pf-pi)/L;

                Ts = v_max/a_max + (e_tot/2)/v_max;
                Tt = (L*a_max+v_max^2)/(a_max*v_max) + e_tot/v_max;

                bound_1 = L-(v_max^2/(2*a_max)) - e_tot/2;
                if bound_1 > L
                    bound_1 = L;
                end

                bound_2 = (v_max^2/(2*a_max)) + e_tot/2;
                if bound_2 < 0
                    bound_2 = 0;
                end

                if bound_1 < bound_2
                    bound_1 = L/2;
                    bound_2 = L/2;
                end
                  
                if (norm(ee_position_h-x_d) >= bound_1)
                    x_dot_d_h = v * a_max*t;
                    x_ddot_d_h = v * a_max;
                elseif (norm(ee_position_h-x_d) < bound_1 && norm(ee_position_h-x_d) >= bound_2)
                    x_dot_d_h = v * v_max;
                    x_ddot_d_h = zeros(3,1);
                elseif (norm(ee_position_h-x_d) < bound_2)
                    x_dot_d_h = - v * a_max*(t-Tt);
                    x_ddot_d_h = - v * a_max;
                end

                m1 = length(x_ddot_d_h);

                %x_ddot_d_h = x_ddot_d_h + 0.001*((x_d - ee_position_h) / norm(x_d - ee_position_h));
                
                x_dot_d = [x_dot_d, x_dot_d_h];
                x_ddot_d = [x_ddot_d, x_ddot_d_h];
                
                
                % SNS solution
                q_ddot_new = SNS_acceleration_multitask(ndof, {m1}, {J1_h}, {J1_dot_h}, {x_ddot_d_h}, bounds, q_h, q_dot_h, T, false);
                q_dot_new = q_dot_h + q_ddot_new*T;
                q_new = q_h + q_dot_h*T + 0.5*q_ddot_new*T^2;

                saturated_acc = false;
                for i=1:ndof
                    if round(q_ddot_new(i),4)>=round(bounds_max_acceleration(i),4) || round(q_ddot_new(i),4)<=round(bounds_min_acceleration(i),4)
                        saturation_counter=saturation_counter+1;
                        saturated_acc = true;
                        break;
                    end
                end
                    
                x_ddot_r = J1_h*q_ddot_new + J1_dot_h*q_dot_new;
                x_dot_r = J1_h*q_dot_new;
                e_ddot = norm(x_ddot_d_h-x_ddot_r);
                e_dot = norm(x_dot_d_h-x_dot_r);
                e_h = e_dot*T+0.5*e_ddot*T^2;

                if saturated_acc
                    e_tot = e_tot + e_h;
                end

                fprintf('==============================================================================\n')

                fprintf('k = %d, norm(ee_position_h - x_d) = %f\n\n', k, norm(ee_position_h-x_d));

                %fprintf('v_max         = ');disp(v_max);
                %fprintf('a_max         = ');disp(a_max);
                %fprintf('pi            = ');disp(pi');
                %fprintf('pf            = ');disp(pf');
                %fprintf('v             = ');disp(v');
                fprintf('t             = ');disp(t);
                fprintf('Ts            = ');disp(Ts);
                fprintf('Tt            = ');disp(Tt);
                fprintf('bound_1       = ');disp(bound_1);
                fprintf('bound_2       = ');disp(bound_2);
                fprintf('ee_position   = ');disp(self.robot.get_ee_position(q_new)');
                fprintf('e_ddot        = ');disp(e_ddot);
                fprintf('e_dot         = ');disp(e_dot);  
                fprintf('e_h           = ');disp(e_h);
                fprintf('e_tot         = ');disp(e_tot);
                fprintf('x_d           = ');disp(x_d');
                fprintf('x_dot_d_h     = ');disp(x_dot_d_h');
                fprintf('x_ddot_d_h    = ');disp(x_ddot_d_h');
                
                fprintf('q_ddot_new    = ');disp(q_ddot_new')
                fprintf('q_dot_new     = ');disp(q_dot_new')
                fprintf('q_new         = ');disp(q_new')

                fprintf('sat_counter   = ');disp(saturation_counter);

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
                ee_position_hm1 = ee_position_h;
        
                % update current cunfiguration to new configuration
                q_h = q_new;
                q_dot_h = q_dot_new;
                q_ddot_h = q_ddot_new;

                J1_h = self.robot.get_J(q_h);
                J1_dot_h = self.robot.get_J_dot(q_h, q_dot_h);
                ee_position_h = self.robot.get_ee_position(q_h);

                if  ((norm(x_dot_r) < self.epsilon || t>=Tt) && t>(Tt-Ts)) || norm(ee_position_h - x_d) < 0.0001
                    fprintf('norm(x_dot_d_h) < self.epsilon           = ');disp(norm(x_dot_d_h) < self.epsilon);
                    fprintf('t>=Tt                                    = ');disp(t>=Tt);
                    fprintf('norm(ee_position_h - x_d) < self.epsilon = ');disp(norm(ee_position_h - x_d) < 0.0001);

                    k = k+1;
                    if k>size(self.path,2)
                        break;
                    end
                    
                    x_d = self.path(1:3,k);
                    pi = ee_position_h;
                    pf = x_d;
                    t = 0;
                    e_tot = 0;
                    saturation_counter = 0;

                    versors{end+1} = [pi, (pf-pi)/norm(pf-pi)];
                end

                t = t+T;

                % Directional error
                e_d = acos(dot(((x_d-ee_position_h)/norm(x_d-ee_position_h)), ((J1_h*q_dot_h)/norm(J1_h*q_dot_h))));

                % save new configuration
                qs = [qs, q_h];
                eds = [eds, e_d];
            end
        end
    end


end