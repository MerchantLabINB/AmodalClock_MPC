function set_axes_limits_f2(ax)
%% Set figures axes and limits
%Panel A
%Panel B
view(ax(3),[-90 0])
ylim(ax(3),[-120 20])
zlim(ax(3),[-30 100])
set(ax(3),'ytick',[-100 -50 0 50 ])
set(ax(3),'ztick',[0 50 100 ])
%Panel C
view(ax(4),[741 6])
xlim(ax(4),[-50 150])
ylim(ax(4),[-80 150])
% zlim(ax(4),[-80 150])
%Panel D
view(ax(5),[50 7.74])
xlim(ax(5),[-50 150])
ylim(ax(5),[-50 120])
zlim(ax(5),[-80 150])

set(ax(5),'xtick',[-100 -50 0 50 100])
set(ax(5),'ytick',[-50 0 50 ])
set(ax(5),'ztick',[-50 0 50 100 150 200 250])

%Panel E
view(ax(1),[-100 -20])%Panel E
xlim(ax(1),[-100 200])
ylim(ax(1),[-100 150])
zlim(ax(1),[-50 150])

set(ax(1),'Ztick',[-50  0 50 100 150])
set(ax(1),'ytick',[-100 -50 0 50 100])
set(ax(1),'xtick',[-100  0 100 200])

%Panel F
view(ax(6),[58,-29])
xlim(ax(6),[0 100])
ylim(ax(6),[-20 80])
zlim(ax(6),[-20 100])

set(ax(6),'xtick',[0 50 100])
set(ax(6),'ytick',[0 40 80])
set(ax(6),'ztick',[-20 0 50 100])

%set labels, ticks and lindwidth
set(ax,'Linewidth',2)
set(ax,'TickDir','out');

xlabel(ax,'PC 1','FontSize',12)
ylabel(ax,'PC 2','FontSize',12)
zlabel(ax,'PC 3','FontSize',12)

