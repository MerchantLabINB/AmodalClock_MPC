function ApplyParams(gcf,MeanColor)
 MeanColor = flip(MeanColor);
    lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
    for i=1:4
        set(lines(i),'Color',MeanColor(i,:),'linewidth',2)
    end
    h = findobj(gca,'Tag','Box');
    
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),MeanColor(j,:),'EdgeColor',MeanColor(j,:),'FaceAlpha',.4);
    end
    
     lines = findobj(gcf, 'type', 'line', 'Tag', 'Box');
    for i=1:4
        set(lines(i),'Color',MeanColor(i,:),'linewidth',1.5)
    end
    
axis(gca,'square')
set(gca,'linewidth',1.5)
% set(gca,'FontSize',12)
set(gca,'TickDir','out');
 set(gca,'box','off');
   