classdef KUKA_LBR_IV
    properties
        % robot model from Robotics System Toolbox
        robot

        % number of degrees of freedom of robot
        ndof
        
        % symbolic variables vector
        q
        q_dot
        q_ddot
        
        % simbolic jacobian matrices
        J
        J_dot

        % symbolic end effector position
        ee_position

        % bounds of robot on position, velocity and acceleration
        bounds_position
        bounds_velocity
        bounds_acceleration

        % inital configuraiton
        q_0
        q_dot_0
        q_ddot_0
        ee_position_0

        % name of robot
        name = 'KUKA_LBR_IV';
    end

    methods
        %% constructur
        function self = KUKA_LBR_IV(q_0, q_dot_0, q_ddot_0)

            % import robot from robotics system toolbox
            self.robot = NaN;
            self.ndof = 7;

            % define symbolic variables
            syms t 
            syms q_1(t) q_2(t) q_3(t) q_4(t) q_5(t) q_6(t) q_7(t)
            syms q_dot_1(t) q_dot_2(t) q_dot_3(t) q_dot_4(t) q_dot_5(t) q_dot_6(t) q_dot_7(t)
            syms q_ddot_1(t) q_ddot_2(t) q_ddot_3(t) q_ddot_4(t) q_ddot_5(t) q_ddot_6(t) q_ddot_7(t)
        
            q_1=q_1(t); q_2=q_2(t); q_3=q_3(t); q_4=q_4(t); q_5=q_5(t); q_6=q_6(t); q_7=q_7(t);
            self.q = [q_1;q_2;q_3;q_4;q_5;q_6;q_7];
            
            q_dot_1=q_dot_1(t); q_dot_2=q_dot_2(t); q_dot_3=q_dot_3(t); q_dot_4=q_dot_4(t); q_dot_5=q_dot_5(t); q_dot_6=q_dot_6(t); q_dot_7=q_dot_7(t);
            self.q_dot = [q_dot_1;q_dot_2;q_dot_3;q_dot_4;q_dot_5;q_dot_6;q_dot_7];

            q_ddot_1=q_ddot_1(t); q_ddot_2=q_ddot_2(t); q_ddot_3=q_ddot_3(t); q_ddot_4=q_ddot_4(t); q_ddot_5=q_ddot_5(t); q_ddot_6=q_ddot_6(t); q_ddot_7=q_ddot_7(t);
            self.q_ddot = [q_ddot_1;q_ddot_2;q_ddot_3;q_ddot_4;q_ddot_5;q_ddot_6;q_ddot_7];

            % define DH matrix alpha, d, a, theta
            DH_matrix = [[ pi/2 0     0 q_1];
                         [-pi/2 0     0 q_2];
                         [-pi/2 0.400 0 q_3];
                         [ pi/2 0     0 q_4];
                         [ pi/2 0.390 0 q_5];
                         [-pi/2 0     0 q_6];
                         [    0 0     0 q_7];];

            % symbolic end effector position: ee_position
            T = {1,self.ndof};
            for i=1:self.ndof
                alpha_i = DH_matrix(i,1);
                d_i = DH_matrix(i,2);
                a_i = DH_matrix(i,3);
                theta_i = DH_matrix(i,4);
                
                T_i = [[cos(theta_i) -sin(theta_i)*cos(alpha_i)  sin(theta_i)*sin(alpha_i)  a_i*cos(theta_i)];
                      [sin(theta_i)  cos(theta_i)*cos(alpha_i)  -cos(theta_i)*sin(alpha_i)  a_i*sin(theta_i)];
                      [           0               sin(alpha_i)                cos(alpha_i)               d_i];
                      [           0                          0                          0                  1];];
            
                T{i} = T_i;
            end
            T_0_7 = T{1}*T{2}*T{3}*T{4}*T{5}*T{6}*T{7};
            self.ee_position = T_0_7(1:3, 4);

            % symbolic robot jacobian: J
            self.J = simplify(jacobian(self.ee_position, self.q));

            % symbolic derivate of robot jacobian: J_dot
            self.J_dot = subs(simplify(diff(self.J, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4), diff(q_5), diff(q_6), diff(q_7)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'),str2sym('q_dot_3(t)'),str2sym('q_dot_4(t)'),str2sym('q_dot_5(t)'),str2sym('q_dot_6(t)'),str2sym('q_dot_7(t)')});
            
            % define bounds
            bounds_max_position = [deg2rad(170), deg2rad(120), deg2rad(170), deg2rad(120), deg2rad(170), deg2rad(120), deg2rad(170)];
            bounds_min_position = -bounds_max_position;
            self.bounds_position = double([bounds_min_position;bounds_max_position]);
            
            bounds_max_velocity = [deg2rad(100), deg2rad(110), deg2rad(100), deg2rad(130), deg2rad(130), deg2rad(180), deg2rad(180)];
            bounds_min_velocity = -bounds_max_velocity;
            self.bounds_velocity = double([bounds_min_velocity;bounds_max_velocity]);
            
            bounds_max_acceleration = [deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300)];
            bounds_min_acceleration = -bounds_max_acceleration;  
            self.bounds_acceleration = double([bounds_min_acceleration;bounds_max_acceleration]);

            % define initial configuration
            self.q_0 = double(q_0);
            self.q_dot_0 = double(q_dot_0);
            self.q_ddot_0 = double(q_ddot_0);
            self.ee_position_0 = double(subs(self.ee_position, self.q, self.q_0));
        end

        %% get Jacobian: J
        function J = get_J(self, q)            
            J = double(subs(self.J, self.q, q));
        end

        %% get derivative of Jacobian: J_dot
        function J_dot = get_J_dot(self, q, q_dot)
            J_dot = double(subs(self.J_dot, [self.q, self.q_dot], [q, q_dot]));
        end

        %% get end effector position: ee_position
        function ee_position = get_ee_position(self, q)
            ee_position = double(subs(self.ee_position, self.q, q));
        end
    end
end