function ax =  PlotKinematics(DistProye,Margangle,StartingPoint,f1)


cmap= ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...   
    1 0 0]);

cols = 5;
rows = 8;
SP = 42:42:168; % #bins for each duration x condition
set(0,'CurrentFigure',f1)
%% Amplitude
ax(1) = subplot(rows,cols,2);
nCond = 4;

for i = nCond:-1:1
    line(1:max(SP),DistProye(:,i), ...
        'Color',cmap(i,:),'LineWidth',1);
end

ymin = round(min(DistProye(:))-5);
ymax = round(max(DistProye(:))+7);
ylim([ymin,ymax])%Articulo  Point 1
line([SP; SP],[ymin ymax],'Color',[1 1 1],'LineStyle','--','LineWidth',1)

%% Angle proyection
ax(2) = subplot(rows,cols,3);

for i = nCond:-1:1
    line(1:168,Margangle(:,i), ...
        'Color',cmap(i,:),'LineWidth',1);
end

ylims = [pi/4 pi];
line([SP; SP],[ylims(1),ylims(2)],'Color',[1 1 1],'LineStyle','--','LineWidth',1)
ylim([ylims(1),ylims(2)])
ax(2) = gca;
%% Position
ax(3) = subplot(rows,cols,4);
for i = 1:4
        hold on;
        line(1:168,StartingPoint(:,i), ...
        'Color',cmap(i,:),'LineWidth',1);
end

ylims = [-70 60];
line([SP; SP],[ylims(1),ylims(2)],'Color',[1 1 1],'LineStyle','--','LineWidth',1)
