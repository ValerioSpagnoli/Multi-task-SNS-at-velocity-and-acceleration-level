classdef KUKA_LBR_IV
    properties

        % number of degrees of freedom of robot
        ndof
        t
        
        % symbolic variables vector
        q
        q_dot
        q_ddot
        
        ee_position

        % simbolic jacobian matrices for end effector
        J_ee
        J_dot_ee

        % symbolic elbow position
        elbow_position

        % simbolic jacobian matrices for elbow
        J_eb
        J_dot_eb

        % bounds of robot on position, velocity and acceleration
        bounds_position
        bounds_velocity
        bounds_acceleration

        % inital configuraiton
        q_0
        q_dot_0
        q_ddot_0
        ee_position_0
        elbow_position_0

        % name of robot
        name = 'KUKA_LBR_IV';
    end

    methods
        %% constructur
        function self = KUKA_LBR_IV(q_0, q_dot_0, q_ddot_0)

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

            DH_matrix = [[ pi/2 0.3105 0 q_1];
                         [-pi/2 0      0 q_2];
                         [-pi/2 0.4000 0 q_3];
                         [ pi/2 0      0 q_4];
                         [ pi/2 0.3900 0 q_5];
                         [-pi/2 0      0 q_6];
                         [    0 0.0780 0 q_7];];



            % symbolic end effector position: ee_position
            T = {1,self.ndof};
            for i=1:self.ndof
                alpha_i = DH_matrix(i,1);
                d_i = DH_matrix(i,2);
                a_i = DH_matrix(i,3);
                theta_i = DH_matrix(i,4);
                
                T_i = [[cos(theta_i) -sin(theta_i)*cos(alpha_i)   sin(theta_i)*sin(alpha_i)  a_i*cos(theta_i)];
                       [sin(theta_i)  cos(theta_i)*cos(alpha_i)  -cos(theta_i)*sin(alpha_i)  a_i*sin(theta_i)];
                       [           0               sin(alpha_i)                cos(alpha_i)               d_i];
                       [           0                          0                          0                  1];];
            
                T{i} = T_i;
            end

            %% end effector position
            T_0_7 = T{1}*T{2}*T{3}*T{4}*T{5}*T{6}*T{7};
            self.ee_position = T_0_7(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_ee = simplify(jacobian(self.ee_position, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_ee = subs(simplify(diff(self.J_ee, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4), diff(q_5), diff(q_6), diff(q_7)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'),str2sym('q_dot_3(t)'),str2sym('q_dot_4(t)'),str2sym('q_dot_5(t)'),str2sym('q_dot_6(t)'),str2sym('q_dot_7(t)')});

            %% elbow position
            T_0_4 = T{1}*T{2}*T{3}*T{4};
            self.elbow_position = T_0_4(1:3, 4);
            
            % symbolic robot jacobian for elbow: J_eb
            self.J_eb = simplify(jacobian(self.ee_position, self.q));

            % symbolic derivate of robot jacobian for elbow: J_dot_eb
            self.J_dot_eb = subs(simplify(diff(self.J_eb, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4), diff(q_5), diff(q_6), diff(q_7)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'),str2sym('q_dot_3(t)'),str2sym('q_dot_4(t)'),str2sym('q_dot_5(t)'),str2sym('q_dot_6(t)'),str2sym('q_dot_7(t)')});
            
            %% define bounds
            bounds_max_position = [deg2rad(170), deg2rad(120), deg2rad(170), deg2rad(120), deg2rad(170), deg2rad(120), deg2rad(170)];
            bounds_min_position = -bounds_max_position;
            self.bounds_position = double([bounds_min_position;bounds_max_position]);
            
            bounds_max_velocity = [deg2rad(100), deg2rad(110), deg2rad(100), deg2rad(130), deg2rad(130), deg2rad(180), deg2rad(180)];
            bounds_min_velocity = -bounds_max_velocity;
            self.bounds_velocity = double([bounds_min_velocity;bounds_max_velocity]);
            
            bounds_max_acceleration = [deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300), deg2rad(300)];
            bounds_min_acceleration = -bounds_max_acceleration;  
            self.bounds_acceleration = double([bounds_min_acceleration;bounds_max_acceleration]);

            %% define initial configuration
            self.q_0 = double(q_0);
            self.q_dot_0 = double(q_dot_0);
            self.q_ddot_0 = double(q_ddot_0);
            self.ee_position_0 = double(subs(self.ee_position, self.q, self.q_0));
            self.elbow_position_0 = double(subs(self.elbow_position, self.q, self.q_0));
        end

        %% get Jacobian: J
        function J = get_J_ee(self, q)            
            J = double(subs(self.J_ee, self.q, q));
        end

        function J = get_J_eb(self, q)            
            J = double(subs(self.J_eb, self.q, q));
        end

        %% get derivative of Jacobian: J_dot
        function J_dot = get_J_dot_ee(self, q, q_dot)
            J_dot = double(subs(self.J_dot_ee, [self.q, self.q_dot], [q, q_dot]));
        end

        function J_dot = get_J_dot_eb(self, q, q_dot)
            J_dot = double(subs(self.J_dot_eb, [self.q, self.q_dot], [q, q_dot]));
        end

        %% get end effector position: ee_position
        function ee_position = get_ee_position(self, q)
            ee_position = double(subs(self.ee_position, self.q, q));
        end

        %% get elbow position: elbow_position
        function ee_position = get_elbow_position(self, q)
            ee_position = double(subs(self.elbow_position, self.q, q));
        end
    end
end