
%%ComputedDistanceDiagonalMatrices2023.m
%%Panels A-C Figure 6


%%Hugo Merchant Sept 5/2023

close all
clear
load ('TSDFMatrixForCorrelSept2023');


computeSOdistance = 1;
Ttask = 4;
TotalsBins = 20;
TSO = 10;
close all
colormap('jet');%('autumn');
cmap = colormap;

Bintime = [22.5 42.5 22.5 42.5];

MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';

SOL(1,:) = 'SO 1';
SOL(2,:) = 'SO 2';
SOL(3,:) = 'SO 3';
SOL(4,:) = 'SO 4';

clear TDifference TDistan
clear TDifferenceSO TDistanSO TDAIinitendSO TDAImiddleSO
clear TCellCorrelSO ProbCorrelQuarter
clear sigcellsquarter  CellCorr ProbCorrelQuarter
clear DiagonalIndexSO DistanceIndexSO

v = 1:10;
C = nchoosek(v,2);
clear tran;


plotlocation = [1 2 3 4 6 7 8 11 12 16];

C = [1	1
    1	2
    1	3
    1	4
    2	2
    2	3
    2	4
    3   3
    3   4
    4   4];

NumPaiwisecomaprisons = numel(C(:,1));

if computeSOdistance
    for task = 2:2%Ttask
        tcels = numel(TBinnedSDF{task,1}(:,1));
        for c = 1:NumPaiwisecomaprisons
            
            tso1 = TBinnedSDF{task,C(c,1)};
            tso2 = TBinnedSDF{task,C(c,2)};
            
            [maxtso1,idxmax] = max(tso1');
            [maxtso1,idxmax2] = sort(idxmax);
            tso1 = tso1(idxmax2,:);
            tso2 = tso2(idxmax2,:);
            
            tit = ['Binned Activity ' [MIL(task,:)] ' Sorted by SO1'];
            f1 = figure(0+task);
            set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
            subplot(4,4,plotlocation(c));
            tsot2 = flip(tso2);
            imagesc(tsot2);
            set(gca, 'XTickLabel',[]);
            set(gca, 'YTickLabel',[]);
            zmax = max(max(tso2))*.8;
            caxis([0 zmax])
            axis square; %square
            
            
            for bin2 = 1:TotalsBins
                for bin = 1:TotalsBins
                    Distan11(bin2,bin) = norm(tso1(:,bin)- tso1(:,bin2));
                end
            end
            
            
            for bin2 = 1:TotalsBins
                for bin = 1:TotalsBins
                    DistanT(bin2,bin) = norm(tso1(:,bin)- tso2(:,bin2));
                end
            end
            
            [Min11,inxMin11] = min(Distan11);
            MinDegrees11 = (inxMin11.*360)./TotalsBins; %Results in deegrees
            MinDegrees11(MinDegrees11> 359.999999) = 360;
            
            [Min12,inxMin12] = min(DistanT);
            MinDegrees12 = (inxMin12.*360)./TotalsBins; %Results in deegrees
            MinDegrees12(MinDegrees12> 359.999999) = 360;
            
            DifferenceT = MinDegrees12 - MinDegrees11;
            DifferenceT(DifferenceT< -180.999999) = DifferenceT(DifferenceT< -180.999999) +360;
            
            
            DiagonalIndexSO(C(c,1),C(c,2)) = {DifferenceT};
            DistanceIndexSO(C(c,1),C(c,2)) = {Min12};
            
            TDifferenceSO{task}(C(c,1),C(c,2)) = sum(abs(DifferenceT));
            TDistanSO{task}(C(c,1),C(c,2)) = mean(Min12);
            TDAIinitendSO(C(c,1),C(c,2)) = round(sum(abs(DifferenceT([1 2 3 4 17 18 19 20]))),0);
            TDAImiddleSO(C(c,1),C(c,2)) = round(sum(abs(DifferenceT(5:16))),0);
            TDistaninitendSO(C(c,1),C(c,2)) = round(mean((Min12([1 2 3 4 17 18 19 20]))),0);
            TDistanmiddleSO(C(c,1),C(c,2)) = round(mean((Min12(5:16))),0);
            
            DAI = round(sum(abs(DifferenceT)),0);
            DAIinitend = round(sum(abs(DifferenceT([1 2 3 4 17 18 19 20]))),0);
            DAImiddle = round(sum(abs(DifferenceT(5:16))),0);
            DI = round(mean(Min12),1);
            
            tit = ['Distance between SOs ' [MIL(task,:)]];
            f1 = figure(5+task);
            set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
            subplot(4,4,plotlocation(c));
            hold on
            imagesc(DistanT);
            axis([1 20 1 20])
            caxis([0 1100])
            plot(1:20,inxMin12,'o', 'MarkerEdgeColor',[0 0 0],'LineWidth', 1,'MarkerFaceColor',[1 1 1],'MarkerSize', 6);
            
            set(gca, 'XTickLabel',[]);
            set(gca, 'YTickLabel',[]);
            
            axis square;
            hold off
            
            CellCorr =  zeros(tcels,tcels);
            
            for cel2 = 1:tcels
                for cel = 1:tcels
                    CellCorr(cel2,cel) = corr(tso1(cel,:)', tso2(cel2,:)');
                end
            end
            
            InxSigCellCorr =  CellCorr>0.5;
            numquartercells = round(tcels/4)-1;
            Tquarter = 4;
            TotalProbNeuQuarter = numquartercells*numquartercells;
            
            
            ncel2 = 1;
            for quar2 = 1:Tquarter
                ncel1 = 1;
                for quar1 = 1:Tquarter
                    sigcellsquarter(quar1,quar2) = sum(sum(InxSigCellCorr(ncel1:ncel1+numquartercells-1, ncel2:ncel2+numquartercells-1)));
                    ProbCorrelQuarter(quar1,quar2) =  sigcellsquarter(quar1,quar2) / TotalProbNeuQuarter;
                    ncel1=ncel1+numquartercells;
                end
                ncel2=ncel2+numquartercells;
            end
            
            TCellCorrelSO{task,C(c,1),C(c,2)} = ProbCorrelQuarter;
            
            tit = ['Probablity of Significant Correlations Between Quarters ' [MIL(task,:)]];
            f1 = figure(10+task);
            set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
            subplot(4,4,plotlocation(c));
            hold on
            imagesc(ProbCorrelQuarter);
            axis([0.5 4.5 0.5 4.5])
            caxis([0 0.4])
            set(gca, 'XTickLabel',[]);
            set(gca, 'YTickLabel',[]);
            axis square;
            hold off
            
            
        end
    end
    
    for task = 2:2%Ttask
        F2 = figure(0+task);
        colormap('jet');
        print(['DischageHeatmap' [MIL(task,:)]],'-dpdf','-painters','-bestfit',  '-r400');%
        saveas(F2,'DischageHeatmapTasksMeanSO.png')
        saveas(F2,'DischageHeatmapTasksMeanSO.fig')
        
        F7 = figure(5+task);
        colormap('jet');
        print(['DistanceSO' [MIL(task,:)]],'-dpdf','-painters','-bestfit',  '-r400');%
        saveas(F7,'DistanceSO.png')
        saveas(F7,'DistanceSO.fig')
        
        F12 =figure(10+task);
        colormap('jet');
        print(['ProbCorrelSO' [MIL(task,:)]],'-dpdf','-painters','-bestfit',  '-r400');%
    end
end



tcels = numel(TBinnedSDF{1}(:,1));
tbins = numel(TBinnedSDF{1}(1,:));
clear t1 MeanSObinned
for task = 1:Ttask
    for so = 1:4%TSO
        t1(so,1:tcels,1:tbins) = TBinnedSDF{task,so};
    end
    tem = mean(t1,1);
    tem2 = reshape(tem, tcels,tbins);
    
    [maxtso1,idxmax] = max(tem2');
    MeanSObinned{task} =  tem2;
end

plotlocation = [1 2 3 4 6 7 8 11 12 16];

C = [1	1
    1	2
    1	3
    1	4
    2	2
    2	3
    2	4
    3   3
    3   4
    4   4];



NumPaiwisecomaprisons = numel(C(:,1));
colormap('jet');%('autumn');
cmap = colormap;
clear TaskDifference TaskDistance
clear sigcellsquarter  CellCorr ProbCorrelQuarter
clear TCellCorreltask TaskDAIinitend TaskTDAImiddle
clear DiagonalIndexTask DistanceIndexTask

for task = 4:4%4:Ttask %just for taggin the figures
    colormap('jet');
    for c = 1:NumPaiwisecomaprisons
        
        
        tso1 = MeanSObinned{C(c,1)};
        tso2 = MeanSObinned{C(c,2)};
        
        
        [maxtso1,idxmax] = max(tso1');
        [maxtso1,idxmax2] = sort(idxmax);
        tso1 = tso1(idxmax2,:);
        tso2 = tso2(idxmax2,:);
        
        tit = ['Mean Binned Activity Across SO Sorted by Diagonal'];
        f1 = figure(15+task);
        set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
        subplot(4,4,plotlocation(c));
        tsot2 = flip(tso2);
        imagesc(tsot2);
        set(gca, 'XTickLabel',[]);
        set(gca, 'YTickLabel',[]);
        zmax = max(max(tso2))*0.8;
        caxis([0 zmax])
        %caxis([20 80])
        axis square;
        
        for bin2 = 1:TotalsBins
            for bin = 1:TotalsBins
                Distan11(bin2,bin) = norm(tso1(:,bin)- tso1(:,bin2));
            end
        end
        
        
        for bin2 = 1:TotalsBins
            for bin = 1:TotalsBins
                DistanT(bin2,bin) = norm(tso1(:,bin)- tso2(:,bin2));
            end
        end
        
        [Min11,inxMin11] = min(Distan11);
        MinDegrees11 = (inxMin11.*360)./TotalsBins; %Results in deegrees
        MinDegrees11(MinDegrees11> 359.999999) = 360;
        [Min12,inxMin12] = min(DistanT);
        MinDegrees12 = (inxMin12.*360)./TotalsBins; %Results in deegrees
        MinDegrees12(MinDegrees12> 359.999999) = 360;
        
        DifferenceT = MinDegrees12 - MinDegrees11;
        DifferenceT(DifferenceT< -180.999999) = DifferenceT(DifferenceT< -180.999999) +360;
        TaskDifference(C(c,1),C(c,2)) = sum(abs(DifferenceT));
        TaskDistance(C(c,1),C(c,2)) = mean(Min12);
        TaskDAIinitend(C(c,1),C(c,2)) = round(sum(abs(DifferenceT([1 2 3 4 17 18 19 20]))),0);
        TaskTDAImiddle(C(c,1),C(c,2)) = round(sum(abs(DifferenceT(5:16))),0);
        
        TDistaninitend(C(c,1),C(c,2)) = round(mean((Min12([1 2 3 4 17 18 19 20]))),0);
        TDistanmiddle(C(c,1),C(c,2)) = round(mean((Min12(5:16))),0);
        
        DiagonalIndexTask(C(c,1),C(c,2)) = {DifferenceT};
        DistanceIndexTask(C(c,1),C(c,2)) = {Min12};
        
        DAI = round(sum(abs(DifferenceT)),0);
        DAIinitend = round(sum(abs(DifferenceT([1 2 3 4 17 18 19 20]))),0);
        DAImiddle = round(sum(abs(DifferenceT(5:16))),0);
        DI = round(mean(Min12),2);
        
        for ii = 12:20
            if  inxMin12(ii) < 10
                inxMin12(ii) = 20-inxMin12(ii);
            end
        end
        
        
        tit = ['Distance between Tasks for the mean SO'];
        f1 = figure(20+task);
        set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
        subplot(4,4,plotlocation(c));
        hold on
        imagesc(DistanT);
        axis([1 20 1 20])
        caxis([0 1100])
        
        plot(1:20,inxMin12,'o', 'MarkerEdgeColor',[0 0 0],'LineWidth', 1,'MarkerFaceColor',[1 1 1],'MarkerSize', 6);
        set(gca, 'XTickLabel',[]);
        set(gca, 'YTickLabel',[]);
        
        axis square;
        hold off
        
        CellCorr =  zeros(tcels,tcels);
        for cel2 = 1:tcels
            for cel = 1:tcels
                CellCorr(cel2,cel) = corr(tso1(cel,:)', tso2(cel2,:)');
            end
        end
        
        InxSigCellCorr =  CellCorr>0.5;
        
        numquartercells = round(tcels/4)-1;
        Tquarter = 4;
        TotalProbNeuQuarter = numquartercells*numquartercells;
        
        if c == 1
            tit = ['Correlation matrix between cells'];
            f1 = figure(30);
            set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
            
            hold on
            imagesc(InxSigCellCorr);
            axis([1 908 1 908])
            caxis([0 1])
            set(gca, 'XTickLabel',[100:100:900]);
            set(gca, 'YTickLabel',[100:100:900]);
            xlabel('Cell #')
            ylabel('Cell #')
            
            axis square;
            colormap('jet');
            hold off
            print(['CorrelatioMatrixA450'],'-dpdf','-painters','-bestfit',  '-r400');%
        end
        
        ncel2 = 1;
        for quar2 = 1:Tquarter
            ncel1 = 1;
            for quar1 = 1:Tquarter
                sigcellsquarter(quar1,quar2) = sum(sum(InxSigCellCorr(ncel1:ncel1+numquartercells-1, ncel2:ncel2+numquartercells-1)));
                ProbCorrelQuarter(quar1,quar2) =  sigcellsquarter(quar1,quar2) / TotalProbNeuQuarter;
                ncel1=ncel1+numquartercells;
            end
            ncel2=ncel2+numquartercells;
        end
        
        TCellCorreltask{task,C(c,1),C(c,2)} = ProbCorrelQuarter;
        
        tit = ['Probablity of Significant Correlations Between Quarters Between Tasks uning the mean SO' ];
        f1 = figure(25+task);
        set(f1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1],'color', 'black','Name', tit);
        subplot(4,4,plotlocation(c));
        hold on
        imagesc(ProbCorrelQuarter);
        axis([0.5 4.5 0.5 4.5])
        caxis([0 0.4])
        
        set(gca, 'XTickLabel',[]);
        set(gca, 'YTickLabel',[]);
        axis square;
        hold off
        
    end
end

F17 = figure(15+task);
colormap('jet');
print(['DischageHeatmapTasksMeanDM'],'-dpdf','-painters','-bestfit',  '-r400');%
saveas(F17,'DischageHeatmapTasksMeanDM.png')
saveas(F17,'DischageHeatmapTasksMeanDM.fig')


F22 = figure(20+task);
colormap('jet');
print(['DistanceTasksMeanDM'],'-dpdf','-painters','-bestfit',  '-r400');%
saveas(F17,'DistanceTasksMeanDM.png')
saveas(F17,'DistanceTasksMeanDM.fig')

figure(25+task);
colormap('jet');
print(['ProbCorrelTasksMeanSO'],'-dpdf','-painters','-bestfit',  '-r400');%

