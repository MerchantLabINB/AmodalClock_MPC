function f = Trajectories_analyses(pathTofolder)
%%Ploting neural trajectories and subspaces for  the four task condition.

% task condition:
%                 Auditory-450ms, Auditory-850ms 
%                 Visual - 450 ms, Visual-850ms 
% for simplicity we call short and long durations respectively.
%each condition has a block of 25 trials.
% Monkey 1  -  8 sessions
% Monkey 2 -  14 sessions.

% We analized the first serial order of the Synchronization task (ST).
% subpaces; separatrix, duration, serial order and modality plane.
% Amodal clock

% Input: Binarized Firing rates: 22  and 42 bins for each duration
% respectively.
%Output: we saved diferents varaibles useful for computing next figure (3)
% of the article.

%%% MERCHANT LAB 2023

% Filename1 = 'Trajectories_LS200_TapControl_1019';
Filename1 = 'FiringRate_binarized';
Outputfile1 = 'Manifold.mat';
VideoWarp = 'mov_dwell_warp'; %We need to compute this file
Videospeed = 'VideoSpeedWarp'; %We need to compute this file
handspeed = 'handspeed_get_index';
fullFileName1 = fullfile(pathTofolder,Filename1);
fullFileName2 = fullfile(pathTofolder,Outputfile1);

%load struct file
load(fullFileName1)

interval = [450,850,450,850];
trials = [25 25 25 25];

if~exist('FR','var')
    FR = {[]};
end

outputStruct{1,1}.ContorlAct = 0;
outputStruct{1}.BinMode = 'Variable';
outputStruct{1,1}.ModalityFr = 'Average';
outputStruct{1,1}.interval = interval;
outputStruct{1,1}.trials = trials;
outputStruct{1,1}.threshold = 0;%4Hz
outputStruct{1,1}.FRControl = FR;
outputStruct{1,1}.sd = Cells;
outputStruct{1}.Projection = 'All'; %Projection
outputStruct{1,1}.AvgCoeff = 1; %Use Average Coefficients
outputStruct{1,1}.coeffPCA =  []; %Coeff;
outputStruct{1,1}.normDataCoeff = 1;%1
outputStruct{1}.PlotProjection = 'Taps';%Projection
outputStruct{1,1}.figurepath = '.\Figures\Figure2';
outputStruct{1,1} = get_trajectoriesPopulation(outputStruct{1});

Manifold = Statistics_Updated(outputStruct);
f = Manifold.figure;
Manifold = rmfield(Manifold,'figure');

if exist(fullFileName2, 'file') == 2
    sprintf('Warning: file  exist:\t%s', Outputfile1)
else
    load(fullfile(pathTofolder,VideoWarp))
    load(fullfile(pathTofolder,Videospeed))
    load(fullfile(pathTofolder,handspeed))

    Manifold.mov_dwell_warp = mov_dwell_warp;
    Manifold.VideoSpeedWarp = VideoSpeedWarp;
    Manifold.HandTrajectoriesSpeed =  HandTrajectoriesSpeed ;
    Manifold.get_index_Movement = get_index_Movement;


    save(fullFileName2,'Manifold')
end
