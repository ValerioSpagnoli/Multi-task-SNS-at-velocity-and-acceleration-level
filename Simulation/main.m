function siml = main()
    clear all; clc; close all;

    level = 'acceleration';
    robot_name = 'Planar4R'; % Planar4R, KUKA_LBR_IV, KUKA_LBR_IIWA_7_R800
    q_0 = NaN;
    q_dot_0 = NaN;
    q_ddot_0 = NaN;
    simulation_step = NaN;
    epsilon = NaN;

    if strcmp(robot_name, 'Planar4R')
        siml = Simulation_Planar4R(level, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon);
    elseif strcmp(robot_name, 'KUKA_LBR_IV') || strcmp(robot_name, 'KUKA_LBR_IIWA_7_R800')
        siml = Simulation_KUKA7R(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon);
    end
end