classdef Planar4R
    properties
        % number of degrees of freedom of robot
        ndof
        t

        % symbolic variables vector
        q
        q_dot
        q_ddot
        
        % symbolic end effector position
        link_1
        link_2
        link_3
        ee_position

        % simbolic jacobian matrices for end effector
        J_link_1
        J_dot_link_1

        J_link_2
        J_dot_link_2

        J_link_3
        J_dot_link_3

        J_ee
        J_dot_ee

        % bounds of robot on position, velocity and acceleration
        bounds_position
        bounds_velocity
        bounds_acceleration

        % inital configuraiton
        q_0
        q_dot_0
        q_ddot_0

        link_1_0
        link_2_0
        link_3_0
        ee_position_0        

        % name of robot
        name = 'Planar4R';
    end

    methods
        %% constructur
        function self = Planar4R(q_0, q_dot_0, q_ddot_0)
            self.ndof = 4;

            % define symbolic variables
            syms t 
            syms q_1(t) q_2(t) q_3(t) q_4(t)
            syms q_dot_1(t) q_dot_2(t) q_dot_3(t) q_dot_4(t)
            syms q_ddot_1(t) q_ddot_2(t) q_ddot_3(t) q_ddot_4(t)
            
            self.t = t;
        
            q_1=q_1(t); q_2=q_2(t); q_3=q_3(t); q_4=q_4(t);
            self.q = [q_1;q_2;q_3;q_4];
            
            q_dot_1=q_dot_1(t); q_dot_2=q_dot_2(t); q_dot_3=q_dot_3(t); q_dot_4=q_dot_4(t);
            self.q_dot = [q_dot_1;q_dot_2;q_dot_3;q_dot_4];

            q_ddot_1=q_ddot_1(t); q_ddot_2=q_ddot_2(t); q_ddot_3=q_ddot_3(t); q_ddot_4=q_ddot_4(t);
            self.q_ddot = [q_ddot_1;q_ddot_2;q_ddot_3;q_ddot_4];

            % define DH matrix alpha, d, a, theta
            DH_matrix = [[0 0 0.25 q_1];
                         [0 0 0.25 q_2];
                         [0 0 0.25 q_3];
                         [0 0 0.25 q_4];];

            
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

            %% link_1 position
            T_0_1 = T{1};
            self.link_1 = T_0_1(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_link_1 = simplify(jacobian(self.link_1, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_link_1 = subs(simplify(diff(self.J_link_1, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'), str2sym('q_dot_3(t)'), str2sym('q_dot_4(t)')});

            %% link_2 position
            T_0_2 = T{1}*T{2};
            self.link_2 = T_0_2(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_link_2 = simplify(jacobian(self.link_2, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_link_2 = subs(simplify(diff(self.J_link_2, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'), str2sym('q_dot_3(t)'), str2sym('q_dot_4(t)')});

            %% link_3 position
            T_0_3 = T{1}*T{2}*T{3};
            self.link_3 = T_0_3(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_link_3 = simplify(jacobian(self.link_3, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_link_3 = subs(simplify(diff(self.J_link_3, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'), str2sym('q_dot_3(t)'), str2sym('q_dot_4(t)')});
            
            %% end effector position
            T_0_4 = T{1}*T{2}*T{3}*T{4};
            self.ee_position = T_0_4(1:3, 4);

            % symbolic robot jacobian for end effector: J_ee
            self.J_ee = simplify(jacobian(self.ee_position, self.q));

            % symbolic derivate of robot jacobian for end effector: J_dot_ee
            self.J_dot_ee = subs(simplify(diff(self.J_ee, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'), str2sym('q_dot_3(t)'), str2sym('q_dot_4(t)')});

            %% define bounds
            bounds_max_position = [2*pi, 2*pi, 2*pi, 2*pi];
            bounds_min_position = -bounds_max_position;
            self.bounds_position = double([bounds_min_position;bounds_max_position]);
            
            bounds_max_velocity = [2, 2, 4, 4];
            bounds_min_velocity = -bounds_max_velocity;            
            self.bounds_velocity = double([bounds_min_velocity;bounds_max_velocity]);
            
            bounds_max_acceleration = [2, 2, 4, 4];
            bounds_min_acceleration = -bounds_max_acceleration;
            self.bounds_acceleration = double([bounds_min_acceleration;bounds_max_acceleration]);

            %% define initial configuration
            self.q_0 = double(q_0);
            self.q_dot_0 = double(q_dot_0);
            self.q_ddot_0 = double(q_ddot_0);

            self.link_1_0 = double(subs(self.link_1, self.q, self.q_0));
            self.link_2_0 = double(subs(self.link_2, self.q, self.q_0));
            self.link_3_0 = double(subs(self.link_3, self.q, self.q_0));            
            self.ee_position_0 = double(subs(self.ee_position, self.q, self.q_0));
        end

        %% get Jacobians link_1: J, J_dot
        function J = get_J_link_1(self, q)            
            J = double(subs(self.J_link_1, self.q, q));
        end
        
        function J_dot = get_J_dot_link_1(self, q, q_dot)
            J_dot = double(subs(self.J_dot_link_1, [self.q, self.q_dot], [q, q_dot]));
        end

        function link_1 = get_link_1(self, q)
            link_1 = double(subs(self.link_1, self.q, q));
        end

        %% get Jacobians link_2: J, J_dot
        function J = get_J_link_2(self, q)            
            J = double(subs(self.J_link_2, self.q, q));
        end
        
        function J_dot = get_J_dot_link_2(self, q, q_dot)
            J_dot = double(subs(self.J_dot_link_2, [self.q, self.q_dot], [q, q_dot]));
        end

        function link_2 = get_link_2(self, q)
            link_2 = double(subs(self.link_2, self.q, q));
        end

        %% get Jacobians link_3: J, J_dot
        function J = get_J_link_3(self, q)            
            J = double(subs(self.J_link_3, self.q, q));
        end
        
        function J_dot = get_J_dot_link_3(self, q, q_dot)
            J_dot = double(subs(self.J_dot_link_3, [self.q, self.q_dot], [q, q_dot]));
        end

        function link_3 = get_link_3(self, q)
            link_3 = double(subs(self.link_3, self.q, q));
        end

        %% getters end effectors: ee_position, J, J_dot
        function J = get_J_ee(self, q)            
            J = double(subs(self.J_ee, self.q, q));
        end

        function J_dot = get_J_dot_ee(self, q, q_dot)
            J_dot = double(subs(self.J_dot_ee, [self.q, self.q_dot], [q, q_dot]));
        end

        function ee_position = get_ee_position(self, q)
            ee_position = double(subs(self.ee_position, self.q, q));
        end
    end
end