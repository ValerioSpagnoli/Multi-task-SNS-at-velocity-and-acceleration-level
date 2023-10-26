function siml = main()
    clear all; clc; close all;

    level = 'acceleration';
    robot_name = 'Planar4R'; % KUKA_LBR_IV, KUKA_LBR_IIWA_7_R800
    q_0 = NaN;
    q_dot_0 = NaN;
    q_ddot_0 = NaN;
    simulation_step = NaN;
    epsilon = NaN;
    round_up = NaN;
    round_point = NaN;

    siml = Simulation(level, robot_name, q_0, q_dot_0, q_ddot_0, simulation_step, epsilon, round_up, round_point);
end