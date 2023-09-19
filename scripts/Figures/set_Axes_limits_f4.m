function set_Axes_limits_f4(ax,Tcells)

    xlabel(ax(1),'Time (3.4 s)')
    ylabel(ax(1),'#Cells')
    set(ax(1))
    set(ax(1))

    set(ax, 'YTick', ([0:100:Tcells]))
    set(ax, 'YTickLabel', ([0:100:Tcells]))
    set(ax, 'XTick', ([]))
    set(ax, 'TickLength', [.01 .01],'linewidth',1.3);            %  this in proportion of x axis%

    set(ax, 'YAxisLocation', 'left')
    set(ax, 'LineWidth', 1,'Color','w','FontSize',10);
    set(ax, 'Box', 'off');
    set(ax, 'TickDir', 'out');                     %  switch side of axis for tick marks

