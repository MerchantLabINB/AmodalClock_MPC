
function f1 = PlotMovingBumpsJul2023Aud(begin,DUr,DischargeRate,peak,PeakMag,task,Serial_order,TapTimes,PeakQuarterNum,f1)

set(0,'CurrentFigure',f1)

cmap2 = colormap;
colormap('summer');
cmap = colormap;
%quartercolors = [1 20 40 64 30 1];
quartercolors = [1 40 120 240 180 1];
MaxDischargeRate = 30;
MaxColor = 64;
ROW_HEIGHT = 1;
Delta = 0;
TapTimes = round(TapTimes);
x_tick_span = TapTimes;%round(mean(diff(TapTimes)));

pre = 250;
if task == 1 | task == 3
    post = 2500;
elseif task == 2 | task == 4
    post = 4200;
end
post = 4500;

MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';

plot_rows = 5;%1
plot_cols = 4;%2


ax1 = subplot(plot_rows, plot_cols, [task,task+4]);

LocationIndex = 0;
cuatro = -35;
inxcolor = 1;


hold on;

TotalCells = numel(begin);

trial_max = 420;%830;%800;%120;%TotalCells;
beg_x = TapTimes;%-starttap;
end_x = TapTimes;%-starttap;
top_y = (0);
bot_y = (-trial_max);
y = [top_y bot_y];
x = [beg_x; end_x];


line(x,y,'Color','w','LineWidth',1);

begin = begin - Delta;
peak = peak - Delta;

colorstep = 64/TotalCells;
tinxcolor = -colorstep+1;



for cellidx = 1:TotalCells
    
    
    LocationIndex = LocationIndex+1;
    tinxcolor = tinxcolor+colorstep;
    inxcolor = round(tinxcolor+colorstep);
    
    if inxcolor > 64
        inxcolor = 64;
    end
    
    quarterC = GetQuarterColor(PeakQuarterNum,cellidx);
    
    
    colorNow = cmap(quartercolors(quarterC),:);%quartercolors(quarterC);
    
    
    if (DischargeRate(cellidx) >MaxDischargeRate)
        DischargeRate(cellidx) = MaxDischargeRate;
    end
    % first trial is top row, build raster by adding trials below it
    row_y_loc = -((LocationIndex-1) * ROW_HEIGHT);
    
    % define y values (rasters drawn top down).
    top_y = (row_y_loc);
    bot_y = (row_y_loc);
    y = [top_y bot_y];
    
    Begt = begin(cellidx);
    Endt = Begt+DUr(cellidx);%round(Activations.ActBeg{duration,seq}(ac) + TeoDur(duration,seq) + Activations.ActDur{duration,seq}(ac));
    Peak = peak(cellidx);%round(Activations.Peak{duration,seq}(ac) + TeoDur(duration,seq));
    if PeakMag(cellidx)>MaxDischargeRate
        PeakMag(cellidx)=MaxDischargeRate;
    end
    
    PeakMagColor = round((PeakMag(cellidx)*MaxColor)/MaxDischargeRate);
    
    x = [Begt Endt];% beg_x end_x];
    xp = Peak;%[peak_x];
    yp = top_y;
    
    
    line(x,y,'Color',colorNow,'LineWidth',0.5);
    plot(xp,yp, 'o', 'MarkerEdgeColor','none','MarkerFaceColor',[cmap2(52,:)],'MarkerSize', 2);
    
    
    if LocationIndex == 1
        axis([pre, post, -trial_max * ROW_HEIGHT,0]);  %  reset ranges for x and y axes
        tTapTimes = TapTimes - TapTimes(1);
        set(gca, 'XTick', TapTimes);        %  reset tick spacing
        set(gca, 'XColor', 'w');                          %  turn off y ticks
        set(gca, 'YColor', 'w');                          %  turn off y ticks
        set(gca, 'TickDir', 'out');                     %  switch side of axis for tick marks
        set(gca, 'TickLength', [.01 .01]);            %  this in proportion of x axis%
        set(gca, 'XTickLabel', tTapTimes, 'FontSize',12,'Color','k');                     %  turn off x-axis tick labeling
        set(gca, 'LineWidth', 1.5');
        set(gca, 'YTick', [-trial_max:50:0]);
        set(gca, 'YTickLabel', [0:50:trial_max]);
        
        if task == 1
            set(gca, 'YTick', [-trial_max:50:0]);                          %  turn off y ticks
        else
        end
        
        title(ax1,[ MIL(task,:) ],'Color','k','FontSize',12);%,'FontWeight','bold')
        
    end
    
    
end

hold off;


