function plotPieQuarterNum(AllPeakQuarterNum)
figure(3)

colormap('summer');
cmap = colormap;
quartercolors = [1 20 40 64 30 1];
labels = {'Initial','Second','Third','Final'};
tt = sum(AllPeakQuarterNum)
tt(5) = [];
pie(tt,labels)
print('PieQuarterNum','-dpdf','-painters','-bestfit',  '-r400');% -fillpage')