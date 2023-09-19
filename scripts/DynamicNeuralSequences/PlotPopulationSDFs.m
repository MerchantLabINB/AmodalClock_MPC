
function [OutSDFM,OutSDSD,OutNcells] = PlotPopulationSDFs(SSDFCell,sdinit,sdend,task,Serial_order,ColorCellGroups,BehTimesSDF,PeakQuarterNum,Tquarter)


colormap('summer');%colormap('jet');%('autumn');
cmap = colormap;
%quartercolors = [1 20 40 64 30 1];
quartercolors = [1 40 120 240 180 1];
MaxDischargeRate = 30;
MaxColor = 64;
ROW_HEIGHT = 1;


row_y_loc = [-100 -200 -300 -400 -500 -600];

TapTimes = BehTimesSDF(task).Taps;
x_tick_span = round(mean(diff(TapTimes)));

pre = -100;
post = 3400;

MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';

plot_rows = 2;
plot_cols = 2;

scrsz = get(0,'ScreenSize');
tit = [MIL(task,:) '_SerialO_'  num2str(Serial_order)  ];%'_Dur_'   num2str(Duration) '_AlignDura_'   num2str(Duration2) ];
f1 = figure(92);
set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'white','Name', tit);

whitebg([0 0 0])

ax1 = subplot(plot_rows, plot_cols, task);

sdinit = sdinit-50;
sdend = sdend-50;


hold on;


beg_x = TapTimes;%-starttap;
end_x = TapTimes;%-starttap;
top_y = (0);
bot_y = (row_y_loc(6));
y = [top_y bot_y];
x = [beg_x; end_x];

line(x,y,'Color',[1 1 1], 'LineWidth',1);
%line(x,y,'Color',[0 0 0], 'LineWidth',1);

psum = 0;
tPeakQuarterNum(1)  = 0;
for q = 1:5
    psum = PeakQuarterNum(q)+psum;
    tPeakQuarterNum(q+1) = psum;
end
LastC = tPeakQuarterNum(end);

Ct = tabulate(ColorCellGroups);
if numel(Ct(:,1)) == 3
    Color2_3 = Ct(2:3,2);
else
    Color2_3 = [0 0];
end


tPeakQuarterNum(6) = tPeakQuarterNum(5) + Color2_3(1);
tPeakQuarterNum(7) = LastC;



for quarter = 1:Tquarter
    
    nnsdf =  tPeakQuarterNum(quarter+1);
    if nnsdf == 0
        SDFq = zeros(1,numel(SSDFCell(1,:)));
        SDFsd = zeros(1,numel(SSDFCell(1,:)));
    elseif nnsdf == 1
        SDFq = (SSDFCell(tPeakQuarterNum(quarter)+1:tPeakQuarterNum(quarter+1),:));
        SDFsd = (SSDFCell(tPeakQuarterNum(quarter)+1:tPeakQuarterNum(quarter+1),:));
    else
        
        SDFq = nanmean(SSDFCell(tPeakQuarterNum(quarter)+1:tPeakQuarterNum(quarter+1),:));
        SDFsd = nanstd(SSDFCell(tPeakQuarterNum(quarter)+1:tPeakQuarterNum(quarter+1),:));
    end
    
    quarterC = GetQuarterColor2(tPeakQuarterNum,tPeakQuarterNum(quarter));
    
    colorNow = cmap(quartercolors(quarterC),:);%quartercolors(quarterC);
    
    if quarterC == 5 %High peak variability between SO
        colorNow = [1 0 0 ];%%red
    end
    
    if quarterC == 6 %%one peak for all sos
        colorNow = [1 1 1];%WHITE
        %colorNow = [0 0 0];
    end
    
    
    x = [sdinit:sdend];
    %     for nc = tPeakQuarterNum(quarter)+1:tPeakQuarterNum(quarter+1)
    %         y = [row_y_loc(quarter)+SSDFCell(nc,:)];
    %          line(x,y,'Color',colorNow,'LineWidth',1.5);
    %     end
    % define y values (rasters drawn top down).
    
    y = [row_y_loc(quarter)+(SDFq*5)];
    
    
    
    
    line(x,y,'Color',colorNow,'LineWidth',1.5);
    
    
    
    if quarter == 1
        axis([pre, post, row_y_loc(6),0]);  %  reset ranges for x and y axes
        set(gca, 'XTick', TapTimes,'FontSize',8);        %  reset tick spacing
        set(gca, 'XColor', 'k');                          %  turn off y ticks
        set(gca, 'YColor', 'k');                          %  turn off y ticks
        set(gca, 'TickDir', 'out');                     %  switch side of axis for tick marks
        set(gca, 'TickLength', [.01 .01]);            %  this in proportion of x axis%
        set(gca, 'XTickLabel', [TapTimes]);                     %  turn off x-axis tick labeling
        set(gca, 'LineWidth', 1.5');
        %   set(gca, 'YTickLabel', [-row_y_loc]);
        %   set(gca, 'YTickLabel', [ ]);
        set(gca, 'YTick', [-row_y_loc]);
        %         set(gca, 'xlabel' (ax2,'Bin')
        if task == 1
            
        else
            %             set(gca, 'YTick', []);
            %             ylabel( ' ');
        end
        
        title(ax1,[ MIL(task,:) ],'Color','k','FontSize',12);%,'FontWeight','bold')
        
    end
    if isempty(SDFq)
    else
        OutSDFM(quarter,:) = SDFq;
        OutSDSD(quarter,:) = SDFsd;
        OutNcells(quarter,:) = nnsdf;
    end
end

hold off;

%print('SDFs4tasks','-dpdf','-painters','-bestfit',  '-r400');% -fillpage')

