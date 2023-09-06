clear all; clc; close all;

%% Define a 4R robot

syms t 
syms q_1(t) q_2(t) q_3(t) q_4(t)
syms q_dot_1(t) q_dot_2(t) q_dot_3(t) q_dot_4(t)
syms q_ddot_1(t) q_ddot_2(t) q_ddot_3(t) q_ddot_4(t)

n = 4;

q_1=q_1(t); q_2=q_2(t); q_3=q_3(t); q_4=q_4(t);
q = [q_1;q_2;q_3;q_4;];

q_dot_1=q_dot_1(t); q_dot_2=q_dot_2(t); q_dot_3=q_dot_3(t); q_dot_4=q_dot_4(t);
q_dot = [q_dot_1;q_dot_2;q_dot_3;q_dot_4;];

q_ddot_1=q_ddot_1(t); q_ddot_2=q_ddot_2(t); q_ddot_3=q_ddot_3(t); q_ddot_4=q_ddot_4(t);
q_ddot = [q_ddot_1;q_ddot_2;q_ddot_3;q_ddot_4;];


bounds_max_position = [2*pi, 2*pi, 2*pi, 2*pi];
bounds_min_position = -bounds_max_position;
bounds_position = double([bounds_min_position;bounds_max_position]);

bounds_max_velocity = [2, 2, 4, 4];
bounds_min_velocity = -bounds_max_velocity;
bounds_velocity = double([bounds_min_velocity;bounds_max_velocity]);

bounds_max_acceleration = [2, 2, 4, 4];
bounds_min_acceleration = -bounds_max_acceleration;  
bounds_acceleration = double([bounds_min_acceleration;bounds_max_acceleration]);

bounds = {bounds_position, bounds_velocity, bounds_acceleration};

simulation_step = 1;

%% Define the robot initial state and the tasks

q_0 = [pi/2; -pi/2; pi/2; -pi/2;];
q_dot_0 = zeros(4,1);
q_ddot_0 = zeros(4,1);

% position in the plane
task_1 = [cos(q_1)+cos(q_1+q_2)+cos(q_1+q_2+q_3)+cos(q_1+q_2+q_3+q_4);
          sin(q_1)+sin(q_1+q_2)+sin(q_1+q_2+q_3)+sin(q_1+q_2+q_3+q_4);];

% moving along y direction for the tip of the second link
task_2 = [sin(q_1)+sin(q_1+q_2)];


J_1_sym = simplify(jacobian(task_1, q));
J_2_sym = simplify(jacobian(task_2, q));

J_dot_1_sym = subs(simplify(diff(J_1_sym, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'),str2sym('q_dot_3(t)'),str2sym('q_dot_4(t)')});
J_dot_2_sym = subs(simplify(diff(J_2_sym, t)), {diff(q_1), diff(q_2), diff(q_3), diff(q_4)}, {str2sym('q_dot_1(t)'), str2sym('q_dot_2(t)'),str2sym('q_dot_3(t)'),str2sym('q_dot_4(t)')});

J_1 = double(subs(J_1_sym, q, q_0));
J_2 = double(subs(J_2_sym, q, q_0));

J_dot_1 = double(subs(J_dot_1_sym, [q, q_dot], [q_0, q_dot_0]));
J_dot_2 = double(subs(J_dot_2_sym, [q, q_dot], [q_0, q_dot_0]));

x_ddot_1 = [-3; -1.5]*10e3;
x_ddot_2 = [1];

%% Pseudoinversion solution (minimum norm solution)

q_dot_MNS = pinv(J_1)*x_ddot_1;
fprintf('Pseudoinversion solution (minimum norm solution): q_dot_MNS = ');disp(q_dot_MNS')

%% Apply SNS multitask algorithm

q_ddot_SNS = SNS_acceleration_multitask(n, {length(x_ddot_1),length(x_ddot_2)}, {J_1, J_2}, {J_dot_1, J_dot_2}, {x_ddot_1, x_ddot_2}, bounds, q_0, q_dot_0, simulation_step, true);
fprintf('Saturation in Null Space solution:                q_ddot_SNS = ');disp(round(q_ddot_SNS', 3));
