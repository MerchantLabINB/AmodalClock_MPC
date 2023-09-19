function Kinematics = ComputeKinematicsParameters(pathToFolder)

Hyperplane = 'Hyperplane_PSO';
ManifoldFile = 'Manifold';

%Temporal Files
load(fullfile(pathToFolder,ManifoldFile))
load(fullfile(pathToFolder,Hyperplane))


Manifold = Create_Struct_Mean(Manifold);


% for point = list%1:size(XYZ,1) 
for point = 1%
   
x = XYZ(point,1);
y = XYZ(point,2);
z = XYZ(point,3);


Manifold.DurPry =[x,y,z]; %Point from hyperplane
Manifold.PointLine = 1; %Point from the duration plane
Manifold.PointMod =  5; %Point from the modality plane
Manifold.mov_dwell = Manifold.mov_dwell_warp;

Kinematics = GetKinematics_metrics(Manifold,'DurPry');
end
