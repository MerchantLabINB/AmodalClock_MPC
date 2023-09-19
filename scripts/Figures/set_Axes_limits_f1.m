function set_Axes_limits_f1(ax)
%Panel A,B,D
%set figure axes
set(ax(1:4),'XLIm',[400 900]);
set(ax(1:4),'XTick',[450 850]);
axis(ax(1:4),'square');
box(ax(3),'off');

xlabel(ax(1:3),'Target interval (ms)')
ylabel(ax(1),'Constant error (ms)')
ylabel(ax(2),'Variability (ms)')
ylabel(ax(4),'AutoCorrelation Lag-1')

ylim(ax(1),[-50 50])
ylim(ax(2),[0 60])
ylim(ax(4),[-.5 .25])

set(ax(1),'YTick',[-50 0 50]);
set(ax(2),'YTick',0:20:60);
set(ax(4),'YTick',-.5:0.25:.25);

%Panel F
xlabel(ax(5),'Time (ms)')
ylabel(ax(5),'Speed (a.u.)')
set(ax(5),'Ytick',0:200:800)
box(ax(5),'off');

%Panel G-H
ylabel(ax(6),'Movement duration (ms)')
ylabel(ax(7),'Dwell duration (ms) ')
set([ax(6),ax(7)],'Ylim',[0 600])
set([ax(6),ax(7)],'YTick',0:200:600)
%set labels, ticks and lindwidth
set(ax,'Linewidth',2)
set(ax,'TickDir','out');