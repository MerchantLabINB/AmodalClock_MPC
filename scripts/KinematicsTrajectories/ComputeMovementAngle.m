function AUC = ComputeMovementAngle(params,Manifold)


% sp = [22,42,22,42];
sp = [42,42,42,42];

% addShift = [28;28;32;29];%Distances
addShift = [32;32;32;31];%Angle

mov_dwell = Manifold.mov_dwell.trace;
binsTaps = Manifold.mov_dwell.mov_star;

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
    binss = binsTaps{cModality} - (addShift(cModality)-1);

    len = sum((~isnan(params(:,cModality))));
    xTime  = [1 : 168]*20;
    set(0,'CurrentFigure',f1)
    ax(cModality)= subplot(2, 2, cModality);
    hold on
    line(xTime,params(:,cModality),'Color',cmap(cModality,:),'LineWidth',1.2);
    
    
    ylims = get(gca,'YLim');
    line([SP(cModality,:); SP(cModality,:)],[0 pi],'Color',[1 1 1],'LineStyle','--')

    hold on
    b =1/3;
    plot([1:len ]*20,(mov_dwell{1,cModality}(addShift(cModality):len + (addShift(cModality)-1))*b))

    xlabel(ax(cModality),'Time (ms)')
    
[portionArea_movs{cModality},movs_segment{cModality} ]= getAUC_movs(Manifold.AngleTrials,binss,cModality,sp);

end

AUC.portionArea_movs = portionArea_movs;
AUC.movs_segment = movs_segment;

set(ax,'xlim',[0 xTime(end)])
set(ax,'XTick',[0:1000:3000])
set(ax,'YLim',[0 pi],'YTick',[0 pi],'YTickLabel',{'0' ,'\pi'})
set(ax,'TickDir','out');
ylabel(ax(1),{'Direction', '(Angle vector) (rad)'})
close