function PlotParamsBehaviorLag(bParam,params)
ErrorStd =params.ErrorStd;
CM = params.CM;
MS1 = params.MS1;
ax = params.ax;
f3 = params.f3 ;
for i = 1:4
    %Behavior
    mPositions = bParam(:,i);
    bParams{i,1} = mPositions;
    hold on;
    SEM(1,i) = std(bParams{i,1})/sqrt(length(bParams{i,1}));
    
    
    [ph1, po]= boundedline(ax,params.intervals{1}(i),mean(mPositions),SEM(1,i),['o'],'cmap',CM(i,:),'alpha','transparency', 0.4);
    errorbar(params.intervals{1}(i),mean(mPositions),SEM(1,i)*ErrorStd,'color',CM(i,:),'linewidth',3);
    
    set(ph1,'MarkerSize',MS1,'MarkerEdgeColor',CM(i,:),'MarkerFaceColor',CM(i,:));
end