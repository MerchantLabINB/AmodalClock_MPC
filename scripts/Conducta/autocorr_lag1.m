function  [LagOutput,trialAuto,AutoCorrelationCond] = autocorr_lag1(Behave)
f1 = figure;
C={Behave(:).ID}';
strTask = {'A450','A850','V450','V850'};
Ind = cell2mat(arrayfun(@(x) find(contains(C,strTask{x})), 1:4,'un',0));
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);    
trialAuto = [];
for Cond = 1:4
Taps = {Behave(Ind(:,Cond)).Movs};

for i=1: length(Ind)
 [autocor,~]= arrayfun(@(x) autocorr(diff(Taps{i}(x,1:5))),1:25,'un',0);
 meanAuto(i,:) = mean(cell2mat(autocor'));
   
 catCorr = cell2mat(autocor');
 trialAuto = [trialAuto; catCorr(:,2)];
 AutoCorr(i,:) =  meanAuto(i,2);
end
% AutoCorrCond(Cond,:) = mean(AutoCorr);
AutoCorrCond(:,Cond) = AutoCorr;
AutoCorrelationCond{Cond} = meanAuto;
meanAuto = []; AutoCorr = [];

set(0,'CurrentFigure',f1)
ax(Cond) = subplot(2,2,Cond);
LData = AutoCorrelationCond{Cond};
imagesc(ax(Cond),LData-mean(LData(:)))
colormap(ax(Cond),bluewhitered(64))
cb2 = colorbar(ax(Cond));%,
title(ax(Cond), strTask{Cond},'FontSize',12)
axis square

end
set(ax([2:4]),'Ytick',[],'Xtick',[])% set(gca,'ztick',-150:50:150)2
set(ax(1),'Ytick',1:4:22)% set(gca,'ztick',-150:50:150)2
set(ax(1),'Xtick',1:4)% set(gca,'ztick',-150:50:150)2

ylabel(ax(1),'# Sessions','FontSize',12)
xlabel(ax(1),' Lag ','FontSize',12)

LagOutput = AutoCorrCond;%reshape(AutoCorrCond,2,2)';
close

% xticks(ax2,[1 2 3 4 5])
% xticklabels(ax2,{'1','2','3','4','5'})
% xticks(ax1,[1 2 3 4 4.5])
% xticklabels(ax1,{'1'})

