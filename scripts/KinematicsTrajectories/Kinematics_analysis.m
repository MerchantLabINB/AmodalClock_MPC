function f1 = Kinematics_analysis(pathToFolder)

%%Ploting Kinematic of neural trajectoriess for  the four task condition.

% task condition:
%                 Auditory-450ms, Auditory-850ms 
%                 Visual - 450 ms, Visual-850ms 
% for simplicity we call short and long durations respectively.
%each condition has a block of 25 trials.
% Monkey 1  -  8 sessions
% Monkey 2 -  14 sessions.


% We analized the dynamic of each prodcued interval:
% We used a gemetric and kinematic approches.
%The kinematics of neural trajectories was characterized using
% the amplitude, angle, and relative position.
% 

% Input: neural trajectories and subspaces.

%%% MERCHANT LAB 2023

FileName1 = 'movements_dwells_sessions';
FileName2 = 'Manifold';

%Temporal Files
% load(fullfile(pathToFolder,FileName5))
load(fullfile(pathToFolder,FileName1))
load(fullfile(pathToFolder,FileName2))


cols = 8;
rows = 5;

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
whitebg('k')

set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);

Kinematics = ComputeKinematicsParameters(pathToFolder);


Amplitude = Kinematics.Amplitude;
Angle =  Kinematics.Angle;
Position = Kinematics.Position;

%Panel A,B,C
ax = PlotKinematics(Amplitude,Angle,Position,f1);
%Panel D,E,F,
% figure;
ax= Kinematics_metrics(Kinematics,f1,ax);
%% Panel H,I,J
trajStruct2 = Create_Mean_Speed_Struct(trajStruct);
trajStruct2(1).avgSpeed =0;
SpeedParams = Get_trace_mean_Speed(Manifold.trajectoriesSO,trajStruct2,1);
ax= mov_dwell_trace_panels(SpeedParams,f1,ax);

%%Panel G 
%Compute Speed of the neural trajectories
% Input:
%       Traj.Trajectories_SO - projection of the neural trajectories
%                        onto the first 3 pc.
SpeedTraj = cellfun(@Trajectory_Speed, Manifold.trajectoriesSO, 'UniformOutput', 0);
ax = Speed_trace_Modalities(SpeedTraj,Manifold.VideoSpeedWarp,Manifold.mov_dwell_warp.trace,f1,ax);

%Panel K - AMSI
Kinematics.trajectories = Manifold.trajectories;
ax = AMSI(Kinematics,f1,ax);
%% Panel L
ax = Correlation_lags(Manifold,f1,ax);
%%
set_Axes_limits_f3(ax)
%%
%%%%%%%%%%
%Para esta funcion esta pendiente obtener VideoSpeedWarp function
%y mov_dwell_warp automaticamente.
% Manifold.VideoSpeedWarp = Traj.VideoSpeedWarp;
% Manifold.mov_dwell_warp = Traj.mov_dwell_warp;

%load('./delete_data/hanspeed_getindex')
% Manifold.HandTrajectoriesSpeed = HandTrajectoriesSpeed;
% Manifold.get_index_Movement = get_index_Movement;
