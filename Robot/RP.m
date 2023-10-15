classdef RP
    properties
        % number of degrees of freedom of robot
        ndof
        t

        % symbolic variables vector
        q
        q_dot
        q_ddot
        
        % symbolic end effector position
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

        % name of robot
        name = 'KUKA_LBR_IIWA_7_R800';
    end

    methods
        %% constructur
        function self = RP(q_0, q_dot_0, q_ddot_0)
            self.ndof = 2;

            % define symbolic variables
            syms t 
            syms q_1(t) q_2(t)
            syms q_dot_1(t) q_dot_2(t)
            syms q_ddot_1(t) q_ddot_2(t)
            
            self.t = t;
        
            q_1=q_1(t); q_2=q_2(t);
            self.q = [q_1;q_2;];
            
            q_dot_1=q_dot_1(t); q_dot_2=q_dot_2(t);
            self.q_dot = [q_dot_1;q_dot_2;];

            q_ddot_1=q_ddot_1(t); q_ddot_2=q_ddot_2(t);
            self.q_ddot = [q_ddot_1;q_ddot_2;];

            % define DH matrix alpha, d, a, theta
            DH_matrix = [[pi/2 0   0 q_1];
                         [0    q_2 0   0];];

            
            %% homogeneus matrices
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
            T_0_2 = T{1}*T{2};
            self.ee_position = T_0_2(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_ee = simplify(jacobian(self.ee_position, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_ee = subs(simplify(diff(self.J_ee, t)), {diff(q_1), diff(q_2)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)')});

            %% define bounds
            bounds_max_position = [2*pi, 1];
            bounds_min_position = [-2*pi, 0.5];
            self.bounds_position = double([bounds_min_position;bounds_max_position]);
            
            bounds_max_velocity = [2.5, 2.5];
            bounds_min_velocity = [-0,8, -0.8];
            self.bounds_velocity = double([bounds_min_velocity;bounds_max_velocity]);
            
            bounds_max_acceleration = [10, 10];
            bounds_min_acceleration = [-8, -8];
            self.bounds_acceleration = double([bounds_min_acceleration;bounds_max_acceleration]);

            %% define initial configuration
            self.q_0 = double(q_0);
            self.q_dot_0 = double(q_dot_0);
            self.q_ddot_0 = double(q_ddot_0);
            self.ee_position_0 = double(subs(self.ee_position, self.q, self.q_0));
        end

        %% get Jacobian: J
        function J = get_J_ee(self, q)            
            J = double(subs(self.J_ee, self.q, q));
        end

        %% get derivative of Jacobian: J_dot
        function J_dot = get_J_dot_ee(self, q, q_dot)
            J_dot = double(subs(self.J_dot_ee, [self.q, self.q_dot], [q, q_dot]));
        end

        %% get end effector position: ee_position
        function ee_position = get_ee_position(self, q)
            ee_position = double(subs(self.ee_position, self.q, q));
        end
    end
end