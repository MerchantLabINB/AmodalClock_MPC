function ax = Kinematics_metrics(Traj,f1,ax)
cmap= ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);

params.f3 = f1;
%%

%Dwell amplitude
amp_dwells = cell2mat(cellfun(@(x) mean(x),Traj.AUCdistances.portionArea_dwell,'un',0)); %PI
pi_params{1,1} = reshape(amp_dwells,4,4)';

%Movement amplitude
pi_movs = cell2mat(cellfun(@(x) mean(x),Traj.AUCMovs.portionArea_movs,'un',0));
pi_params{1,2} = reshape(pi_movs,4,4)';


% Variability of Trajectories (Position)
TrajPI = Traj.PositionTrials;
% TrajPI = Traj.AUCPositionDwell.stdPosition;
VarPosition = cell2mat(cellfun(@(x) mean(x),TrajPI  ,'un',0)); %PI
pi_params{1,4} = reshape(VarPosition,4,4)';


%% Dwell - Movs
pi_params{1,3}  = ComputeAutoCorrelationTrials(Traj.AUCdistances.portionArea_dwell);
%%
for i=1:4
    [h,p,ci,stats] = ttest(pi_params{1,3}(i,:));
    p_params(i,1) = p;
    tValue(i,1) = stats.tstat;
    h_params(i,1) = h;
end
set(0,'CurrentFigure',f1)

cols = 5;
rows = 8;


%From Behav Array
params.intervals{1} = [450,850,450,850];
params.CM = [0,0.7490,1; 0 0 1; 1 .5020 0;1 0 0];
params.MS1 = 5;
params.ErrorStd = 2;

%index for subplots
idxplots{2} = [7,6,9,10];

for k=2
    for j = 1:length(idxplots{k})
        
        if( j==1)
            params.ErrorStd = 3;
            %     elseif((k == 3 && j==1)|| (k == 3 && j==2) || (k == 3 && j==2))
        elseif(j==2)
            params.ErrorStd = 2;
        elseif(j==3)

            params.ErrorStd = 1;%.5
        else
            params.ErrorStd = 3;
        end
            
            ax( j + 3) = subplot(rows,cols,idxplots{k}(j));
            params.ax = ax( j + 3);
            PlotParamsBehavior(pi_params{1,j},params)

    end
end

hold on
x=1:900;
y=0;
plot(ax(6),x,y*ones(size(x)),'Color','w','LineStyle','--')
%Article
% set(ax(2,1:2),'ylim',[0 1300],'ytick',[0:400:1200])
% set(ax(2,3),'ylim',[-.5 .25],'ytick',[-.5 0 .25])
% set(ax(2,4),'ylim',[15 40],'ytick',[15:10:35])

% set(ax(2,1:2),'ylim',[0 1800],'ytick',[0:400:1600])
% set(ax(2,3),'ylim',[-.5 .25],'ytick',[-.5 0 .25])
% set(ax(2,4),'ylim',[10 25],'ytick',[10,15,20,25])


%%
% Amp_mov = Traj.AUCMovs.portionArea_movs;
% Amp_dwell = Traj.AUCdistances.portionArea_dwell;
% 
% this_stat = KinematicsStatistics(Amp_mov,Amp_dwell);
% 
% %Estadística actual 2023
% this_stat.autocorrelation = KinematicsStatistics({pi_params{1,3}'});
% this_stat.Position = KinematicsStatistics(TrajPI);

% this_stat = BehaviorStaLagstistics(bParams,Lags);
%2023
%[this_stat.Stats,this_stat.PostHoc] = Statistics_repeatedMeasures(bParams,Lags);









