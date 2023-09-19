function f1 = NeuralSequences_analyses(pathToFolder)

%%Ploting heatmaps neural activity of cells in the four task condition with
%%significan MI effect in at least one of the four key parameters of the
%%experiment: duration, modality, seral order and elapsed time

%%% MERCHANT LAB 2023

FileName1 = 'NeuralSequences';

%Temporal Files
load(fullfile(pathToFolder,FileName1))



f1 = figure(1);%('units','normalized','outerposition',[0 0 1 1]);

set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', 'k')

PosterFigurePopALLHM4T2023Fig_2242(OutData,f1);

function  PosterFigurePopALLHM4T2023Fig_2242(OutData,f1)



MinZ = 0;%Min and Max color scale on normalized activity
MaxZ = 85;

TTask = 4;

MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';

NumBins = [88 168 88 168];

ZV450 = OutData.tDZscoreV450;%V450 D
ZV850 = OutData.tDZscoreV850;%V850 D
ZA450 = OutData.tDZscoreA450;%A450 M
ZA850 = OutData.tDZscoreA850;%A850 M

for task = 1:TTask

    eval_str = [' tResponse =  Z' MIL(task,:) ';'];
    eval(eval_str)

    [tem1,sx] = max(tResponse');
    eval_str = [' [IdxT,indexConditions' MIL(task,:) '] = sort(sx);'];
    eval(eval_str)

end

set(0,'CurrentFigure',f1)
hold on
for task = 1:TTask

    eval_str = [' tResponse =  Z' MIL(task,:) ';'];
    eval(eval_str)

    tSelResp= tResponse;

    [tem1,sx] = max(tSelResp');
    [IdxT,indexConditions] = sort(sx);

    indexConditions = indexConditions';
    ttSelResp= tSelResp(indexConditions,:);
    tbin = numel(ttSelResp(1,:));
    ttSelResp(:,2:tbin+1) = ttSelResp;
    Tcells = numel(indexConditions);

    %%%Short visual selective
    ax(task) = subplot(1,4,task);
    Fline = NumBins(task)/4;
    vbin= linspace(Fline,NumBins(task),4)+.2;
    imagesc(ax(task),ttSelResp)
    line(ax(task),[vbin;vbin],[0,Tcells],'Color','w','linewidth',0.5, 'lineStyle','-')
    title(ax(task),[ MIL(task,:) ],'Color','w','FontSize',12,'FontWeight','bold')

    [cmin(1), cmax(1)] = caxis;
    AxOne   = get(gca,'Parent');
    colormap(AxOne ,'jet')

    ax1.CLim = [MinZ,MaxZ];

end
set_Axes_limits_f4(ax,Tcells)

