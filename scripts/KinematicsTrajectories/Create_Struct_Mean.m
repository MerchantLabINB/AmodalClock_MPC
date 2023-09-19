function Manifold = Create_Struct_Mean(Manifold)
numPc = 3;
lenBin =42;

trajarray = Manifold.trajectories;
trajarray = cellfun(@(x) x(1:numPc,:),trajarray,'un',0);

traj_resh = cellfun(@InterpTrajectory, trajarray , 'UniformOutput', 0);
TrajSO =cellfun(@(x) permute(reshape(x', lenBin, 4, numPc), [1 3 2]),traj_resh, 'un',0);

PryTrajectories = cellfun(@(x) getPryTrialTrajs(Manifold,x), traj_resh,'un',0 );
PryTrajectories_SO = cellfun(@(x) permute(reshape(x', lenBin, 4, numPc), [1 3 2]),PryTrajectories, 'un',0);

%% Cell Array
for i=1:4
temp_Traj = trajarray(:,i)';
auxTemp = nanmean(reshape(cell2mat(temp_Traj), [size(temp_Traj{1,1},1), size(temp_Traj{1,1},2), size(temp_Traj,2)]), 3);
temp_traj{1} = auxTemp; 

if(i==1 || i==3)
  MeanTrj = cellfun(@InterpTrajectory, temp_traj(1) , 'UniformOutput', 0);
else
    MeanTrj{1} = auxTemp;
end

% lenBin =42;
MeanTrajNanSO{i} = permute(reshape(MeanTrj{1}', lenBin, 4, numPc), [1 3 2]);%4

end
% [1 3 2 4]
MeanTrajCellArray{1} = MeanTrajNanSO{1};
MeanTrajCellArray{2} = MeanTrajNanSO{2};
MeanTrajCellArray{3} = MeanTrajNanSO{3};
MeanTrajCellArray{4} = MeanTrajNanSO{4};
Manifold.KiNeT_Cells = MeanTrajCellArray;
Manifold.KineT_TrajPry = PryTrajectories_SO;
Manifold.KineT_Traj = TrajSO;

%%
for i=1:4
DurTraj{i} = permute(reshape(Manifold.DurTraj{i}', lenBin, 4, numPc), [1 3 2]);
end

MeanDurTraj{1} = DurTraj{1};
MeanDurTraj{2} = DurTraj{2};
MeanDurTraj{3} = DurTraj{3};
MeanDurTraj{4} = DurTraj{4};
Manifold.KineT_Dur = MeanDurTraj;
