function CartonParams(R_Dur,R_Mod,R_Tim,f1)
ModC(2,:) = [255,69,0]/255;	%orangered
ModC(1,:) = [36, 84, 163]/255;

% DurC(2,:) = [128,128,0]/255; %olive
% DurC(1,:) = [0,128,0]/255 ; %green
DurC(1,:) = [255,255,0]/255;
DurC(2,:) = [218,165,32]/255;
numMark = 8;
numMark2 = 7;

% whitebg('k')

set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);
set(0,'CurrentFigure',f1)
h = plot3(R_Dur(1, 1:25), R_Dur(2, 1:25), R_Dur(3, 1:25), 'd','color',DurC(1,:));
set(h,'MarkerSize',numMark2,'LineStyle', 'none','MarkerFaceColor',DurC(1,:))
hold on

h = plot3(R_Dur(1, 26:50), R_Dur(2, 26:50), R_Dur(3, 26:50), 'd','color',DurC(2,:));
set(h,'MarkerSize',numMark2,'LineStyle', 'none','MarkerFaceColor',DurC(2,:))
h  = plot3(R_Mod(1, 1:25), R_Mod(2, 1:25), R_Mod(3, 1:25), 'Marker','*','color',ModC(1,:) );
set(h,'MarkerSize',numMark,'LineStyle', 'none')
h = plot3(R_Mod(1, 26:end), R_Mod(2, 26:end), R_Mod(3, 26:end),'color',ModC(2,:), 'Marker','*');
set(h,'MarkerSize',numMark,'LineStyle', 'none');


plot3(R_Tim(1,: ), R_Tim(2, :), R_Tim(3, :), 'w','LineWidth',.3)
binTimeTaps = sort([42:42:168*25]);

hold on; 
plot3(R_Tim(1,binTimeTaps),R_Tim(2,binTimeTaps),R_Tim(3,binTimeTaps), '.','markersize', 40, 'color', [ 0 1 0])
axis equal


