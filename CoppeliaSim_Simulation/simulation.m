function simulation()
    disp('sim')
    clear_all();
    [sim, simID]= start_simulation();
    if simID > -1
        main(sim, simID);
    else
        disp('Failed connecting to remote API server.');
    end
    end_simulation(sim, simID);
end


function clear_all()
    clear all; clc; close all; return;
end
    

function [sim, simID] = start_simulation()
    disp('Program started.');
    sim=remApi('remoteApi');
    sim.simxFinish(-1); 
    simID=sim.simxStart('127.0.0.1',19999,true,true,5000,5);
end


function end_simulation(sim, simID)
    sim.simxGetPingTime(simID);        
    sim.simxFinish(simID);
    disp('Program ended.');
end


function main(sim, simID)
    n = 7; % number of joints (KUKA LBR4 has 7 joints)
    joints = zeros(1,n);

    for i=1:n
        [res, joint] = sim.simxGetObjectHandle(simID,sprintf('LBR4p_joint%d', i),sim.simx_opmode_oneshot_wait);
        if res ~= sim.simx_return_ok 
            fprintf('Error: cannot handle "simxGetObjectHandle" remote API for LBR4p_joint1.\n');
            return;   
        end
        joints(i) = joint;
    end

    pause(2);

    while(sim.simxGetConnectionId(simID) ~= -1)
        %t = sim.simxGetLastCmdTime(simID)/1000.0;
        %if (t > 150) 
        %break; 
        %end

        for i=1:length(joints)
            [joint_position, joint_velocity, joint_acceleration] = getJointState(sim, simID, joints(i));

            fprintf('joint_%d position: %f [rad], %f [deg]\n', i, joint_position, rad2deg(joint_position));
            fprintf('joint_%d velocity: %f [rad/s], %f [deg/s]\n', i, joint_velocity, convangvel([joint_velocity],'deg/s','rad/s'));
            fprintf('joint_%d acceleration: %f [rad/s^2], %f [deg/s^2]\n\n', i, joint_acceleration, convangacc([joint_acceleration],'deg/s^2','rad/s^2'));
        end
        fprintf('----------------------------------------------------------------------------\n')
    end
end




function [joint_position, joint_velocity, joint_acceleration] = getJointState(sim, simID, joint)
    [res_p, ~,joint_position, ~, ~]=sim.simxCallScriptFunction(simID, '/LBR4p', sim.sim_scripttype_childscript,'getJointPosition',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
    [res_v, ~,joint_velocity, ~, ~]=sim.simxCallScriptFunction(simID, '/LBR4p', sim.sim_scripttype_childscript,'getJointVelocity',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
    [res_a, ~,joint_acceleration, ~, ~]=sim.simxCallScriptFunction(simID, '/LBR4p', sim.sim_scripttype_childscript,'getJointAcceleration',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
    if res_p ~= sim.simx_return_ok 
        fprintf('Error: cannot handle "simxCallScriptFunction" with "getJointPosition" script\n.');
        joint_position = -1.0;
    elseif res_v ~= sim.simx_return_ok 
        fprintf('Error: cannot handle "simxCallScriptFunction" with "getJointVelocity" script\n.');
        joint_velocity = -1.0;
    elseif res_a ~= sim.simx_return_ok 
        fprintf('Error: cannot handle "simxCallScriptFunction" with "getJointAcceleration" script\n.');
        joint_acceleration = -1;
    end
end