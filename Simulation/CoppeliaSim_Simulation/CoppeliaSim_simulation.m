classdef CoppeliaSim_simulation
    properties (SetAccess=private)
        robot
        joint_positions
    end

    methods
        function self = CoppeliaSim_simulation(robot, joint_positions)
            self.robot = robot;
            self.joint_positions = joint_positions;

            sim=remApi('remoteApi');
            sim.simxFinish(-1); 
            clientID=sim.simxStart('127.0.0.1',19999,true,true,5000,5);

            if clientID == -1
                disp('Failed connecting to remote API server.');
                sim.simxGetPingTime(clientID);        
                sim.simxFinish(clientID);
            else
                disp('Connection established. Simulation started.')
            end 

            self.main(sim, clientID);

            sim.simxGetPingTime(clientID);        
            sim.simxFinish(clientID);
            disp('Simulation ended. Connection destroyed.')
        end

        %% main function
        function main(self, sim, clientID)
            n = self.robot.ndof;
            robot_name = self.robot.name;
            fprintf('Robot: %s\n', robot_name);
            
            if strcmp(robot_name,'KUKA_LBR_IIWA_7_R800')
                robot_name = 'LBRiiwa7R800';
            elseif strcmp(robot_name,'KUKA_LBR_IV')
                robot_name = 'LBR4p';
            end

            joints = zeros(1,n);
            for i=1:n
                if strcmp(robot_name,'LBRiiwa7R800')
                    joint_name = strcat('LBR_iiwa_7_R800', sprintf('_joint%d',i));
                elseif strcmp(robot_name, 'LBR4p')
                    joint_name = strcat('LBR4p', sprintf('_joint%d',i));
                end

                [res, joint] = sim.simxGetObjectHandle(clientID,joint_name,sim.simx_opmode_oneshot_wait);
                if res ~= sim.simx_return_ok 
                    fprintf('Error: cannot handle "simxGetObjectHandle" remote API for %s.\n',joint_name);
                    return;   
                end

                joints(i) = joint;
            end
        
            pause(2);
        
            i = 1;
            while(sim.simxGetConnectionId(clientID) ~= -1 && i < size(self.joint_positions,2))
                joint_position = self.joint_positions(:,i);
                [~,~,~,~,~]=sim.simxCallScriptFunction(clientID, strcat('/', robot_name), sim.sim_scripttype_childscript,'setConfiguration',joints,joint_position,'',[],sim.simx_opmode_oneshot_wait);
                i = i+1;
            end
        end

        %% get joint state
        function [joint_position, joint_velocity, joint_acceleration] = getJointState(sim, simID, robot_name, joint)
            [res_p, ~,joint_position, ~, ~]=sim.simxCallScriptFunction(simID, strcat('/',robot_name), sim.sim_scripttype_childscript,'getJointPosition',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
            [res_v, ~,joint_velocity, ~, ~]=sim.simxCallScriptFunction(simID, strcat('/',robot_name), sim.sim_scripttype_childscript,'getJointVelocity',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
            [res_a, ~,joint_acceleration, ~, ~]=sim.simxCallScriptFunction(simID, strcat('/',robot_name), sim.sim_scripttype_childscript,'getJointAcceleration',[joint],[],'',[],sim.simx_opmode_oneshot_wait);
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
    
    end
end