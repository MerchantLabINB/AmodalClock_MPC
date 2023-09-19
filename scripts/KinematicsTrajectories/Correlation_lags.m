function ax = Correlation_lags(ManifoldTime,f1,ax)
% function params = Correlation_lags(ManifoldTime)

%% Correlation analysis using Speed trajectories
% and hand movement speed.
SP = repmat([ 7 23 39 55; 11 34 58 80],2,1);

cols = 5;
rows = 8;


Cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);

%Get index of the movement epoch and collapse them
Cond_indx = ManifoldTime.get_index_Movement;
MovementSpeedMean_nan = ManifoldTime.HandTrajectoriesSpeed;
aux_mean = cellfun(@(x) mean(x),MovementSpeedMean_nan(2,1:4),'un',1);
aux_mean_hand  = cellfun(@(x) mean(x),MovementSpeedMean_nan(1,1:4),'un',1);

set(0,'CurrentFigure',f1)
params.plot = true;
idxPlot = [nan,nan,36,39];
for idxNumlag=8
    for Cond = 1:4 %numConditions

        aux  =MovementSpeedMean_nan{2,Cond};% - aux_mean(Cond);
        y0= [nan(1,idxNumlag),aux(1:end-idxNumlag)]';

        Corr_signal1 =ManifoldTime.HandTrajectoriesSpeed{1,Cond}(Cond_indx{Cond})';%Trajectories
        Corr_signal1_norm = Corr_signal1;%-aux_mean_hand(Cond); 
        Corr_signal2 = y0(Cond_indx{Cond});%Trajectories
        Corr_signal2_norm = Corr_signal2 *2 ;%scale it by factor

        if(params.plot && (Cond ==3 || Cond == 4))
            
            ax(14+Cond) = subplot(rows,cols,[idxPlot(Cond):idxPlot(Cond)+1]);
            hold on;
            yyaxis right
            plot(Corr_signal1_norm,'Color',[150 128 0]/255,'LineWidth',1)
            yyaxis left
            plot(Corr_signal2_norm,'Color',Cmap(Cond,:),'LineWidth',1)
            ylim([0 250])
            line([SP(Cond,:); SP(Cond,:)],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--','LineWidth',1)
            
        end
        
        [Corr_cond(idxNumlag,Cond),Pval(idxNumlag,Cond),] = corr(Corr_signal1,Corr_signal2);
        Corr_signal1 = [];
        Corr_signal1 = [];
    end
end
params.correlation = Corr_cond;
params.pval = Pval;

% set(ax,'Tickdir','out')
% set(ax,'box', 'off')
% set(ax,'ylim',[-50 200])


