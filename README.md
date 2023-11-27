# Multi task Saturation in the Null Space (SNS) method at acceleration and velocity level

This repository contain the framework for computing the Saturation in the Null Space method at acceleration and velocity level, using preemptive prioritization strategy.

The reference papers are:
- F. Flacco, A. De Luca and O. Khatib, "Prioritized multi-task motion control of redundant robots under hard joint constraints," 2012 IEEE/RSJ International Conference on Intelligent Robots and Systems, Vilamoura-Algarve, Portugal, 2012, pp. 3970-3977, doi: 10.1109/IROS.2012.6385619. [Link](https://ieeexplore.ieee.org/document/6385619).
- F. Flacco, A. De Luca and O. Khatib, "Motion control of redundant robots under joint constraints: Saturation in the Null Space," 2012 IEEE International Conference on Robotics and Automation, Saint Paul, MN, USA, 2012, pp. 285-292, doi: 10.1109/ICRA.2012.6225376. [Link](https://ieeexplore.ieee.org/document/6225376).

## Description

The SNS method is a Null Space method to perform prioritized multi task motion control of redudant robot under hard joint constraints, using the preemtive prioritization approach.
Proceeds, iteration by iteration, by discarding the commands that exceed the hard bounds of a task with a specific priority, and reintroducing them at their saturated level, by projection in the Null Space of a suitable Jacobian associated to the already considered tasks (higher priority tasks). Moreover, if the command exceed also with the redistributed effort, a task scaling approach is used to keep the joints into their bounds. 

## Simulation

All simulation are performed at acceleration level, using a sampling time of 1 ms.

### Simulation on Planar 4R

- Simulation 1: 
  - task 1: trajectory tracking task of the hexagon
- Simulation 2: 
  - task 1: trajectory tracking task of the hexagon plus (higher priority),
  - task 2: keep link 1 parallel to the x-axis (lower priority).
- Simulation 3: 
  - task 1: trajectory tracking task of the hexagon plus (higher priority),
  - task 2: self motion dumping, using a configuration space task (lower priority).

| Simulation 1                                                        | Simulation 2                                                                | Simulation 3                                                                |                                    
| -------------------------------------------------------------------- | ----------------------------------------------------------------------------| ----------------------------------------------------------------------------|   
|![Alt Text](/Simulation/Results/Planar4R/task1/Video/gif/cropped.gif) | ![Alt Text](/Simulation/Results/Planar4R/task1_task2/Video/gif/cropped.gif) | ![Alt Text](/Simulation/Results/Planar4R/task1_task3/Video/gif/cropped.gif)|


### Simulation on KUKA LWR IV
- Simulation 1: 
  - task 1: trajectory tracking task of the segment.
- Simulation 2: 
  - task 1: trajectory tracking task of the segment (higher priority),
  - task 2: self motion dumping, using a configuration space task (lower priority).

| Simulation 1                                                        | Simulation 2                                                                |
| ------------------------------------------------------------------- | ----------------------------------------------------------------------------|
|![Alt Text](/Simulation/Results/KUKA_LBR_IV/task1/Video/gif/cropped.gif) | ![Alt Text](/Simulation/Results/KUKA_LBR_IV/task1_task2/Video/gif/cropped.gif) |


## Directory tree
```sh
Multi-task-SNS-at-velocity-and-acceleration-level
├── SNS
│   ├── SNS_velocity_multitask.m
│   └── SNS_accleration_multitask.m
├── Simulation
│   ├── Results
│   ├── Simulations
│   │   ├── Simulation_KUKA.m
│   │   └── Simulation_Planar4R.m
│   ├── MATLAB_Simulation.m
│   └── main.m
└── Robot   
    ├── KUKA_LWR4
    │   ├── meshes
    │   ├── KUKA_LBR_IV.urdf
    │   └── KUKA_LBR_IV.m    
    ├── Planar4R    
    │   ├── Planar4R.urdf
    │   └── Planar4R.m    
    └── KUKA_LBR_IIWA_7R.m
```

## Usage

### Clone the repository
```sh
git clone https://github.com/ValerioSpagnoli/Multi-task-SNS-at-velocity-and-acceleration-level.git
```

### Run the simulation
1. Add the repository to the MATLAB path.
2. Open ``main.m`` and set the variables for a costum simulation:
   1. level: <velocity, acceleration>
   2. robot_name = <Planar4R, KUKA_LBR_IV, KUKA_LBR_IIWA_7_R800> 
   3. q_0: initial configuration of the robot
   4. q_dot_0: initial velocity of the robot;
   5. q_ddot_0 = initial acceleration of the robot;
   6. simulation_step: sampling time of the simulation;
   7. epsilon: the tollerance with which are reached the points of the path;

    **Note**: A value can be set ``NaN``. In this case the simulation will be run with the default parameter. The default parameters are configurated into the simulation files ``Simulation_KUKA.m`` and ``Simulation_Planar4R.m``.

3. Run the file ``main.m`` using the following command in the MATALB console:
    ```
    sim = main();
    ```
4. Show the simulation results using the following command:
    ```
    matlab_simulation = MATLAB_simulation(sim);
    ```
    **Note**: that by setting the ``folder_path`` variable into the file ``MATLAB_simulation.m`` you can save all the plots in a personal folder.

