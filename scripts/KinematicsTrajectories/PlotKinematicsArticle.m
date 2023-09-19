function  PlotKinematicsArticle(DistProye,Margangle,StartingPoint,cmap)
f1 = figure;
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);
set(0,'CurrentFigure',f1)
cols = 3;
rows = 3;

ax = subplot(rows,cols,4);
xLabel = ' Phase (rad)';
nCond = 4;
BinSize = [42:42:168];
numPC = 3;
for i = nCond:-1:1
    line(1:max(BinSize),DistProye(:,i), ...
        'Color',cmap(i,:));
end
SP = 42:42:168;%252
% ylims = [25,150]; %Articulo point 1

ymin = round(min(DistProye(:))-5);
ymax = round(max(DistProye(:))+7);

ylim([ymin,ymax])%Articulo  Point 1
% line([SP; SP],get(gca, 'YLim'),'Color',[1 1 1],'LineStyle','--')
line([SP; SP],[ymin ymax],'Color',[1 1 1],'LineStyle','--')

xlabel(ax,xLabel)
ylabel({'Amplitude (a.u.)'})

set(ax,'YTick',[50 100 150])
% set(ax,'YTick',[0 50 100])

set(ax,'XTick',BinSize)
set(ax,'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'})
set(ax,'TickDir','out');
% axis square

%% Angle proyection
ax = subplot(rows,cols,5);

for i = nCond:-1:1
    line(1:168,Margangle(:,i), ...
        'Color',cmap(i,:));
end
SP = 42:42:168;
ylims = [pi/4 pi];
line([SP; SP],[ylims(1),ylims(2)],'Color',[1 1 1],'LineStyle','--')
ylim([ylims(1),ylims(2)])
xlabel(ax,xLabel)
ylabel({'Angle (rad)'})

set(ax,'XTick',BinSize)
set(ax,'YLim',[pi/4 5*pi/6],'YTick',[pi/4  pi/2 5*pi/6],'YTickLabel',{'\pi/4' ,'\pi/2','5\pi/6'})
set(ax,'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'})
set(ax,'TickDir','out');
% axis square

%% Position
 hold on;
ax = subplot(rows,cols,6);
for i = 1:4%nCond:-1:1
        line(1:168,squeeze(StartingPoint(:,i,numPC)), ...
        'Color',cmap(i,:));
end

ylims = [-70 60];
% ylims = [-40 20];

SP = 42:42:168;
line([SP; SP],[ylims(1),ylims(2)],'Color',[1 1 1],'LineStyle','--')
set(ax,'YLim',ylims,'YTick',[-50 0 50]) %[-50 0 50]

xlabel(ax,xLabel)
ylabel('Position (a.u.)')

set(ax,'XTick',BinSize)
set(ax,'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'})
set(ax,'TickDir','out');
% axis square
