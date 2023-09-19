function AUC = compute_dwell(params,ManifoldTime)




sp = [42,42,42,42]; %Phase

addShift = [28 28 32 29];

% addShift = [27;27;32;27];%28 28 32 29
% addShift = [34 29 35 26]; %New data 2023_1019 cells
% 
%28 28 32 29
mov_dwell = ManifoldTime.mov_dwell.trace;
binsTaps = ManifoldTime.mov_dwell.mov_star;
binsdwell = ManifoldTime.mov_dwell.dwell_Star;


f1 = figure;
whitebg('k')
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);


cmap= ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);


SP(1,:) = [42:42:168];
SP(2,:) = [42:42:168];



SP = repmat(SP,2,1) *20;

for cModality = 1: 4
    
    binss = binsdwell{cModality} - (addShift(cModality)-1);

    len = sum((~isnan(params(:,cModality))));
    xTime  = [1 : 168]*20;
    set(0,'CurrentFigure',f1)
    ax(cModality)= subplot(2, 2, cModality);
    hold on
    line(xTime,params(:,cModality),'Color',cmap(cModality,:),'LineWidth',1.2);
    
    
    ylims = get(gca,'YLim');

    hold on
    plot([1:len ]*20,(mov_dwell{1,cModality}(addShift(cModality):len + (addShift(cModality)-1))*50)+0.2)

    ylim([0.2,ylims(2)])
    xlabel(ax(cModality),'Time (ms)')
    new_trace = mov_dwell{1,cModality}(addShift(cModality):len + (addShift(cModality)-1));
   
[portionArea_dwell{cModality},TotalArea_dwell{cModality},PeaktoPeak_Amp{cModality},dwell_segment{cModality}] = getAUC(ManifoldTime.DPryT,binss,cModality,sp);
end

AUC.portionArea_dwell = portionArea_dwell;
AUC. TotalArea_dwell = TotalArea_dwell;
AUC.PeaktoPeak_Amp = PeaktoPeak_Amp;
AUC.dwell_Segment = dwell_segment;

set(ax,'xlim',[0 xTime(end)])
set(ax,'ylim',[0 120])
set(ax,'XTick',[0:1000:3000])
set(ax,'YTick',[0,50,120])
set(ax,'TickDir','out');
ylabel(ax(1),'Amplitude (a.u.)')
close