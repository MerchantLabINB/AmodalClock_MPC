function Kinematics = GetKinematics_metrics(Manifold,Proyection)
% X should be an N (dimensions/neurons/factors) x C (conditions) x T (time
% points) matrix. Conditions should be ordered in increasing duration.
% Conditions for which the length is less than the max length should be
% padded with NaN values. The length of each condition should be the same
% for all N.
%

%Abraham Betancourt Vera 2021
%%%%%%MERCHANT LAB%%%%%%%%%%%%

scoresM = Manifold.KiNeT_Cells;
LinePoint = Manifold.PointLine; 
PointMod = Manifold.PointMod;


if(strcmp('DurPry',Proyection))
FitPoints  = Manifold.DurPry;
scores = Manifold.KineT_Dur;
else
FitPoints = Manifold.FitLine;
scores = Manifold.KiNeT_Cells;
end

PlanePoints = Manifold.ModPoints;

cmap= ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...   
    1 0 0]);


SO =4;
binarrTaps{1} = reshape(1:42 * SO,42,SO);
binarrTaps{2} = reshape(1:42 * SO,42,SO);
binarrTaps{3} = reshape(1:42 * SO,42,SO);
binarrTaps{4} = reshape(1:42* SO,42,SO);

%% Compute Kinematics Parameters (Amplitude, Angle, Position) onto
%Average neural trajectories.
[DistProye,Margangle,StartingPoint,~,~] = KinematicsParamsProyection(LinePoint,PointMod,FitPoints,scores,scoresM,PlanePoints,binarrTaps);

for i =1:4
% Kinematics.Position{i} = squeeze(StartingPoint(:,i,1));
Kinematics.Position(:,i) = squeeze(StartingPoint(:,i,3));

end
%% Compute Kinematics Parameters (Amplitude, Angle, Position) by
% trials over 4 conditions of the task.
sTrials = Manifold.KineT_TrajPry;
sMTrials = Manifold.KineT_Traj;

[Traj.PhaseTaps,Traj.Ampl,Traj.angle,Traj.SP,Traj.PositionP,Traj.AmplTrials] = KinematicsParamsProyectionTrials(LinePoint,PointMod,FitPoints,sTrials,sMTrials,PlanePoints,binarrTaps,DistProye);

% Compute metrics (dwell & movemente amplitude, angle and variability of
% the position) in kinematics parameters by trials.
Manifold.AmplitudeTrials = Traj.Ampl;
Manifold.AngleTrials = Traj.angle;
Manifold.PositionTrials= Traj.PositionP  ;

%Output
Kinematics.AmplitudeTrials = Traj.Ampl;
Kinematics.AngleTrials = Traj.angle;
Kinematics.PositionTrials= Traj.SP;  
Kinematics.AmplitudeTrialsAMSI= Traj.AmplTrials; 

Kinematics.Amplitude = DistProye;
Kinematics.Angle = Margangle;


Kinematics.AUCdistances = computedwellAmplitude(DistProye,Manifold);
Kinematics.AUCMovs = computeMovementAmplitude(DistProye,Manifold);

Kinematics.AUCangleDwell = ComputedwellAngle(Margangle,Manifold);
Kinematics.AUCangleMovs = ComputeMovementAngle(Margangle,Manifold);

