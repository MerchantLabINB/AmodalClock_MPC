function f1 = Behavior_analyses(pathToFolder)
%%Ploting behavioral parameters for  the four task condition.

% task condition:
%                 Auditory-450ms, Auditory-850ms 
%                 Visual - 450 ms, Visual-850ms 
%each condition has a block of 25 trials.
% Monkey 1  -  8 sessions
% Monkey 2 -  14 sessions.
% We analized the first serial order of the Synchronization task (ST).
% behavioral parameters; Constant error, Variability,Autocorrelation
% We analized the monkeys behavior during the ST by measuring
% the tapping movement kinematics

%%% MERCHANT LAB 2023

Filename1 = 'Behave';
Filename2 = 'Pop.mat';
Filename3 = 'mov_dwell_array.mat';


fullFileName1 = fullfile(pathToFolder,Filename1);
fullFileName2 = fullfile(pathToFolder,Filename2);
fullFileName3 = fullfile(pathToFolder,Filename3);


load(fullFileName1)
load(fullFileName2)


if exist(fullFileName3, 'file')
   load(fullFileName3)

else
  % File does not exist.
  sprintf('Warning: file does not exist:\n%s', fullFileName2)
  sprintf('Running:\t%s', 'dwellSpeed_analysis')
  [mov_dwell_array,~,trajStruct] = mov_dwell_Speed_analysisPop(Pop);
  sprintf('Saving arrays:\t%s', 'mov_dwell_array,movements_dwells_sessions')

  save(fullfile(pathToFolder,'mov_dwell_array'),'mov_dwell_array')
  save(fullfile(pathToFolder,'movements_dwells_sessions'),'trajStruct')
  % warningMessage = sprintf('Warning: file does not exist:\n%s', fullFileName2);
  % uiwait(msgbox(warningMessage));
end

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
whitebg('k')
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);
ax = Conducta_Modelo_2021(Behave,f1); %Panel B,C,D
%Figure 1, Panel G,H
ax= mov_dwell_figure(mov_dwell_array,Pop,f1,ax);
set_Axes_limits_f1(ax)