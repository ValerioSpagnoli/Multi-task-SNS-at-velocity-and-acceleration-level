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

J_1_sym = simplify(jacobian(task_1, q));
J_1 = double(subs(J_1_sym, q, q_0));

x_dot_1 = [-4; -1.5];

%% Pseudoinversion solution (minimum norm solution)

q_dot_MNS = pinv(J_1)*x_dot_1;
fprintf('Pseudoinversion solution (minimum norm solution): q_dot_MNS = ');disp(q_dot_MNS')

%% Apply SNS multitask algorithm

q_dot_SNS = SNS_velocity_multitask(n, {length(x_dot_1)}, {J_1}, {x_dot_1}, bounds, q_0, simulation_step, false);
fprintf('Saturation in Null Space solution:                q_dot_SNS = ');disp(round(q_dot_SNS', 3));
