function set_Axes_limits_f3(ax)

BinSize = 42:42:168;
%% axis square
xlabel(ax(1),' Phase (rad)')
ylabel(ax(1),'Amplitude (a.u.)')
ylabel(ax(2),'Angle (rad)')
ylabel(ax(3),'Position (a.u.)')

set(ax(1),'YTick',[50 100 150])
set(ax(2),'YLim',[pi/4 5*pi/6],'YTick',[pi/4  pi/2 5*pi/6],'YTickLabel',{'\pi/4' ,'\pi/2','5\pi/6'})
set(ax(3),'YLim',[-70 60],'YTick',[-50 0 50])

set(ax(1:3),'XTick',BinSize)
set(ax(1),'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'})
set(ax(1:3),'TickDir','out');
%%
set(ax(4:5),'ylim',[0,1250],'YTick',[0:400:1200])
set(ax(6),'ylim',[-.5,.25],'YTick',[-.5 0 .25])
set(ax(7),'ylim',[15, 40],'YTick',[15 25  35])


xlabel(ax(5),'Produced Interval (ms)')
ylabel(ax(4),'Dwell amplitude')
ylabel(ax(5),'Movement amplitude')
ylabel(ax(6),{'Autocorrelation Lag-1', 'Dwell amplitude'})
ylabel(ax(7),'Variability of Position')

%%
set(ax(8:9),'ylim',[0,80])
set(ax(8),'xtick',[22:22:110],'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'});
set(ax(9),'xtick',[42:42:210],'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'});
ax(8).XAxis.Visible = 'off';
xlabel(ax(9),'Time (bin)')
ylabel(ax(8),'Speed (a.u.)')  

%%
set(ax(10),'xlim',[0 40])
set(ax(11),'xlim',[0 160])

set(ax(10),'ylim',[30 70],'ytick',[30 50 70]);%40
set(ax(11),'ylim',[10 70],'ytick',[10 30 50 70]);%40


set([ax(12),ax(13)],'Ylim',[10 65])
set([ax(12),ax(13)],'YTick',15:15:60)


ylabel(ax(10),'Mov Speed (a.u.)')
ylabel(ax(11),'dwell Speed (a.u.)')

ylabel(ax(12),'Speed (a.u.)')
xlabel(ax(12),'Movements')
ylabel(ax(13),'Speed (a.u.) ')
xlabel(ax(13),'dwells')

set(ax(10:13),'Tickdir','out')
set(ax(10:13),'box', 'off')
axis(ax(10:13),'square')


%%  AMSI

set(ax(14:16),'xlim',[0,168])
set(ax(16),'xtick',42:42:168,'XTickLabel',{'2\pi','4\pi','6\pi','8\pi'})

ax(14).XAxis.Visible = 'off';
ax(15).XAxis.Visible = 'off';

%%
set(ax(17),'xlim',[0,64])
set(ax(18),'xlim',[0,94])
set(ax(17:18),'ylim',[0 250])
set(ax(17:18),'yTick',50:100:250)


xlabel(ax(17:18),'Time (bin)')
ylabel(ax(17),'Speed hand (a.u.)')

yyaxis right
ylabel(ax(18),'Speed (a.u.)')
set(ax(17:18), 'LineWidth', 1,'FontSize',10);
set(ax, 'LineWidth', 1,'FontSize',10);






