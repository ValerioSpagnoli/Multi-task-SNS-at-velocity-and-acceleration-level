classdef Simulation
    properties
        path
        robot
        simulation_step
        joint_positions
    end

    methods
        function self = Simulation(robot, simulation_step)
            self.robot = robot;
            self.simulation_step = simulation_step;

            self.path = self.create_path(2);

            self.joint_positions = self.run_simulation();
        end

        %% create path
        function path = create_path(self, n_cycle)
            poly_0 = [[0; 0; 0.2], [0; 0.1732; 0.1], [0; 0.1732; -0.1], [0; 0; -0.2], [0; -0.1732; -0.1], [0; -0.1732; 0.1]];
            center = [0.1; 0.35; 0.8235];
            
            points = [];
            for n=1:n_cycle
                for i=1:6
                    points = [points, [poly_0(:,i)]];
                end
            end
            path = [points, poly_0(:,1)]+center;
        end

        %% run simulation
        function qs = run_simulation(self)

            ndof = self.robot.ndof;
            bounds = {self.robot.bounds_position, self.robot.bounds_velocity, self.robot.bounds_acceleration};
        
            q_0 = self.robot.q_0;
            q_dot_0 = self.robot.q_dot_0;
            q_ddot_0 = self.robot.q_ddot_0;
            ee_position_0 = self.robot.ee_position_0;
            J_0 = self.robot.get_J(q_0);
            J_dot_0 = self.robot.get_J_dot(q_0, q_dot_0);
        
            T = self.simulation_step;
            kp = 10;
            kd = 0.1;
            
            % current configuration (time h)
            q = q_0;
            q_dot = q_dot_0;
            q_ddot = q_ddot_0;
            ee_position = ee_position_0;
            J = J_0;
            
            % configuration at time h-1
            q_prev = q_0;
            q_dot_prev = q_dot_0;
            q_ddot_prev = q_ddot_0;
            ee_position_prev = ee_position_0;
            J_prev = J_0;
            
            % all configurations during task
            qs = [q_0];
        
            k = 1;
            while true
                
                x = round(self.path(1:3,k),4);
                m = {length(x)};
        
                %fprintf('k = %d, norm(ee_position - x) = %f\n', k, norm(ee_position-x));
            
                if norm(ee_position - x) < 0.01
                    k = k+1;
                    x = self.path(k);            
                end
                if k == length(self.path)+1
                    break;
                end
        
                % compute new configuration
                V = round(kp*norm(x-ee_position) - kd*norm(J_prev*q_dot_prev),4);
                x_dot = round(V * ((x-ee_position)/norm(x-ee_position)),4);
                q_dot_new = SNS_velocity_multitask(m, ndof, {J}, {x_dot}, bounds, q, T, false);
                
                % step forward
                dq = q_dot_new*T;
                q_new = q + dq;
        
                % set previous configuration to current configuration ...
                q_prev = q;
                q_dot_prev = q_dot;
                J_prev = J;
        
                % update current cunfiguration to new configuration
                q = q_new;
                q_dot = q_dot_new;
                J = self.robot.get_J(q);
                ee_position = self.robot.get_ee_position(q);
        
                % save new configuration
                qs = [qs, q];
            end
        end
    end


end