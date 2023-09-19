function ax = Speed_trace_Modalities(Speed,MovementSpeed,binMov,f1,ax)

%% Movements Speed Average
Movement_Modality = arrayfun(@(col) vertcat(MovementSpeed{:, col}), 1:size(MovementSpeed, 2), 'un', 0);
MovementSpeedMean =cellfun(@(x) mean(x)',Movement_Modality,'un',0);
%%

cols = 5;
rows = 8;

CM = [0,0.7490,1; 0 0 1; 1 .5020 0;1 0 0]; %dodgerBlueScheme
set(0,'CurrentFigure',f1)

tt = [0:22:110;0:42:210];

TapsBins(1,:) =[1,(22:22:88),109];
TapsBins(2,:) = [1,(42:42:168),209];
TapsBins(3,:) =[1,(22:22:88),109];
TapsBins(4,:) =[1,(42:42:168),209];




%Speed of Trajecties
% Speed = cellfun(@Trajectory_Speed, Trajectory, 'UniformOutput', 0);
for cModality = 1: 4
    hold on;

    Taps = TapsBins(cModality,:);

    Speed_Modality = [Speed{:, cModality}]';
    Speed_Mean = mean(Speed_Modality);
    

    MovementSpeedMean{1,cModality} = MovementSpeedMean{1,cModality}(16:end);
    binMov{1,cModality} = binMov{1,cModality}(16:end);
    
    if(cModality==1 || cModality==3)
       ax(8) = subplot(rows,cols,[11:12]);
        hold on;
            xlim([0 109])%Varaible
    else
         ax(9) = subplot(rows,cols,[16:17]);
                hold on;
    xlim([0 209])%Varaible

    end
    plot(1: length(Speed_Mean), Speed_Mean, 'color',CM(cModality,:), 'linewidth', 3)
    plot([Taps; Taps],[0 80],'--','Color',[1 1 1],'LineWidth',1.2)

end
