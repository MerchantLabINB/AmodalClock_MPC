function Manifold_ = Statistics_Updated(outputStruct)



% Avgtrajarray = outputStruct{1, 1}.AvgTraj;
trajarray = reshape(outputStruct{1, 1}.traj,25,4);
sdf = cellfun(@(x) x', outputStruct{1, 1}.sd,'un',0 );
sdf = reshape(sdf,25,4);


% outputStruct{1, 1}.Bin = [180 340 180 340]; %10ms
if(~isempty(outputStruct{1, 1}.ControltrajFR{1}))
    trajarrayFr = reshape(outputStruct{1, 1}.ControltrajFR,25,4);
    Pre = size(outputStruct{1, 1}.ControltrajFR{1, 1},2);
else
    trajarrayFr = {[]};
    Pre = 10; %25 %0
end

params.FitLine = true;
params.dims = [1 2 3];% C(comb,:);
params.tasks{1}= 'Blocks';
params.intervals{1} = outputStruct{1, 1}.interval;

%-----Create array of movemnts-----%
meanMovs= arrayfun(@(x) 0: outputStruct{1,1}.Bin(x)/4: outputStruct{1, 1}.Bin(x),1:4,'UniformOutput',false);
meanMovs = cellfun(@(x) [1 x(2:end)] , meanMovs, 'un', false);

%Cut trajectories according to the First and Fifth movement.
%Remove projection from pre and post movement and the last interval.

BinS = size(trajarray{1,1},2) - (size(trajarray{1,1},2)-Pre) + 1 ; %110 PI - 10 Post
BinL = size(trajarray{1,1},2) - (size(trajarray{1,1},2)-Pre) + 1 ; %110 PI - 10 Post
BinC = repmat([BinS,BinL],1,2)   ;
%     BinC = [1 1 1 1];
endBin = (BinC -1) + outputStruct{1, 1}.Bin;


params.Phase = 0;



numCond = size(trajarray,2);
meanTraj = {1};
addBin = outputStruct{1, 1}.Bin/outputStruct{1, 1}.SO;
%-----Average of Trajectories-----%
for i=1:numCond
    temp_Traj = trajarray(:,i)';
    temp_TrajSDF = sdf(:,i);
    auxTemp = mean(reshape(cell2mat(temp_Traj), [size(temp_Traj{1,1},1), size(temp_Traj{1,1},2), size(temp_Traj,2)]), 3);
    %   AverageTraj(i) = cellfun(@(x) x(:,BinC:endBin(i)),temp_AvgTraj,'un',0);

    %Uncomment to normal Analysis
    meanTraj(i) = cellfun(@(x) x(:,BinC(i):endBin(i)),{auxTemp},'un',0);
    trajectories(:,i) =  cellfun(@(x) x(:,BinC(i):endBin(i)),temp_Traj,'un',0);
    Precontrol(:,i) =  cellfun(@(x) x(:,1:BinC(i)-1),temp_Traj,'un',0);
    Postcontrol(:,i) =  cellfun(@(x) x(:,endBin(i)+1:endBin(i)+10),temp_Traj,'un',0);
    sdfTrials(:,i) =  cellfun(@(x) x(:,BinC(i):endBin(i)),temp_TrajSDF,'un',0);

    %     trajectoriestoSpeed(:,i) =  cellfun(@(x) x(:,1:endBin(i)+addBin(i)),temp_Traj,'un',0);
    trajectoriestoSpeed(:,i) = temp_Traj;
    trajectoriesSO(:,i) =  cellfun(@(x) x(:,BinC(i):endBin(i)+addBin(i)),temp_Traj,'un',0);
end

%%
%%  MANIFOLD
TapsTrials = TapsReps(trajectories(:),4); %mov 7 bins forward
Manifold.TapsReps = TapsTrials;
Manifold.movs = meanMovs;
Manifold.Meantrajectories = meanTraj;%change by meanTraj2
Manifold.trajectories = trajectories;%Change by trajectories
Manifold.SDF = sdfTrials;
params.Manifold = Manifold;

TapsTrials = TapsReps(trajectories(:),4); %mov 7 bins forward
Manifold.TapsReps = TapsTrials;

%% Trajectories Params
outputStruct{1}.meanTraj =meanTraj;
outputStruct{1}.intervals =outputStruct{1}.interval;
outputStruct{1}.meanMov = meanMovs;
outputStruct{1}.dims=params.dims;
% f = figure;
f = figure('units','normalized','outerposition',[0 0 1 1]);

ax(1) = subplot(2,3,5);
%% Panel E - Serial Order Plane---------%
if(strcmp(outputStruct{1,1}.BinMode,'Variable'))
    Var_SO = SerialOrder_Plane(Manifold,f);
end

%% Variance and Planes
R = cellfun(@InterpTrajectory, Manifold.trajectories, 'UniformOutput', 0);
R2 = cellfun(@InterpTrajectory, Manifold.trajectories, 'UniformOutput', 0);
MA = R2(:,[1,2]);
MV = R2(:,[3,4]);
R_ModA = [MA{:}]';
R_ModV = [MV{:}]';

R = reshape([R{:}], [3, 168, 25, 2, 2]);
R_Tot = reshape(R, [3, 168*25*2*2]);

%% Variance Planes
R_Dur = reshape(mean(R, 5), 3, 168*25*2);
R_Mod = reshape(mean(R, 4), 3, 168*25*2);

figure
[Vp.DurPlanes,VarExp.DurPlanes,VarPlane(1)] = subspacePlane(R_Dur',1,'Duration'); %Duracion
close
Manifold.V_Dur = Vp.DurPlanes ;
%%  Projections for Kinematics Parameters
if(strcmp(outputStruct{1,1}.BinMode,'Variable'))

    [TrajectoriesDur,TapsPry,TapsRespPry] = getProyections(Manifold);
    [TrajectoriesDurV,TapsPryV,TapsRespPryV] = getProyectionsVariable(Manifold);
    Manifold.DurTraj = TrajectoriesDur;
    Manifold.DurPry = TapsPry;
    Manifold.DurTrajTime = TrajectoriesDurV;
    Manifold.DurPryTime = TapsPryV;
end

%% Panel A - Population Trajectories
ax(2) = subplot(2,3,1);
plotRoutes3D_Updated(params,outputStruct{1}.meanTraj, outputStruct{1}.meanMov,f);

%% Panel B - Duration Plane
params.Plane = [ 0 0 0];
DurTaps{1} = [1,42:42:168];
ax(3) = subplot(2,3,2);
plotRoutes3D_Updated(params,TrajectoriesDur, repmat(DurTaps,1,4),f);

%% Panel C - Duration Plane & separatrix %
ax(4) = subplot(2,3,3);
plotRoutes3D_Updated(params,outputStruct{1}.meanTraj, outputStruct{1}.meanMov,f);
subspacePlane(R_Dur',1,'Duration'); %Duracion

if(params.FitLine) %separatrix
    [xfit,yfit,zfit,Vp.Line] = FittingLine3dSpace(TapsTrials(:,1),TapsTrials(:,2),TapsTrials(:,3));
end
Manifold.FitLine = [xfit,yfit,zfit];

%% Panel D - Modality Plane
ax(5) = subplot(2,3,4);
plotRoutes3D_Updated(params,outputStruct{1}.meanTraj, outputStruct{1}.meanMov,f);
figure;
[Vp.ModPlanes,VarExp.ModPlanes,VarPlane(2)] = subspacePlane(R_Mod',1,'ModalityPlane');%Modalidad
close
[~,~,XYZ] = ModalityPlane(gcf,Manifold);



Manifold.ModPoints = XYZ;
Manifold.PreControl = Precontrol;
Manifold.PostControl = Postcontrol;%change by meanTraj2
% save(fullfile(savepath,'Pop'),'Manifold')

%% Statistics Variance
Var_Tot = trace(cov(R_Tot'));
R_Dur = reshape(mean(mean(R, 5), 2), [3, 25*2]);
Var_Dur = trace(cov(R_Dur'));
R_Mod = reshape(mean(mean(R, 4), 2), [3, 25*2]);
Var_Mod = trace(cov(R_Mod'));
R_Tim = reshape(mean(mean(R, 4), 5), [3, 25*168]);
Var_Tim = trace(cov(R_Tim'));

Var_Mix = Var_Tot - (Var_Dur+Var_Mod+Var_Tim);

TimeDur = [Var_Tot - (Var_Tim + Var_Dur)] / Var_Tot ;
TimeMod = [Var_Tot - (Var_Tim + Var_Mod)]/Var_Tot;
DurMod = [Var_Tot -  (Var_Dur + Var_Mod)]/Var_Tot;

%% Explained Variance by Manifold
[Vp.Tot,VarExp.Tot,VarAcum(1)] = subspacePlane(R_Tot',0);
[Vp.Dur,VarExp.Dur,VarAcum(2)] = subspacePlane(R_Dur',0,'Duration'); %Duracion
[Vp.Mod,VarExp.Mod,VarAcum(3)] = subspacePlane(R_Mod',0);%Modalidad
%
[Vp.ModA,VarExp.ModA,VarAcum(4)] = subspacePlane(R_ModA,0);%Modalidad Auditiva
[Vp.ModV,VarExp.ModV,VarAcum(5)] = subspacePlane(R_ModV,0);%Modalidad Visual
[Vp.Time,VarExp.Time,VarAcum(6)] = subspacePlane(R_Tim',0);%Modalidad Visual

VarAcum(7) = VarExp.Dur(1); %SerialOrder

%Angle
Angle(1) =getAnglePlane(Vp.Dur,Vp.Line); %Duration - Line
Angle(2) =getAnglePlane(Vp.ModA,Vp.ModV); %ModA - ModV
Angle(3) =getAnglePlane(Vp.ModA,Vp.Dur);  %ModA - Duration
Angle(4) =getAnglePlane(Vp.ModV,Vp.Dur);  %ModV - Duration
%% Panel F - Amodal Clock
ax(6) = subplot(2,3,6);
CartonParams(R_Dur,R_Mod,R_Tim,f)

Params.Variance = Vp;
Params.VarExp = VarExp;
Params.VarAcum = VarAcum;
Params.angle = Angle;

%% Output File to next figure (3)
Manifold_.trajectories = Manifold.trajectories;
Manifold_.trajectoriesSO = trajectoriesSO;
Manifold_.FitLine = Manifold.FitLine;
Manifold_.ModPoints = Manifold.ModPoints;
Manifold_.DurTraj = Manifold.DurTraj;
Manifold_.DurPry = Manifold.DurPry;
Manifold_.V_Dur = Manifold.V_Dur ;
% save(fullfile(savepath,'Params'),'Params')
%% Trajectories (Speed)
% figParams.savepath = savepath;
% figParams.Speedfig = Speedfig;
% figParams.Vectorfig = Vectorfig;
% % Speed_Modalities(trajectories,figParams)
% % Speed_Modalities(trajectoriesSO,figParams)
% Speed_Modalities(trajectoriestoSpeed,figParams)
% save(fullfile(savepath,'trajectories'),'trajectoriestoSpeed')
% save(fullfile(savepath,'trajectoriesSO'),'trajectoriesSO')
set_axes_limits_f2(ax)
Manifold_.figure = f;