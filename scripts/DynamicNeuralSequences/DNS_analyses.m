function f1 = DNS_analyses(pathTofolder)
%%Figure 5 Panel A
%%ploting and saving all the parameters of activation periods within a neural sequence from spike train poisson analysis
%% for the 4 tasks conditions

%%% MERCHANT LAB 2023

f1 = figure('units','normalized','outerposition',[0 0 1 1]);
whitebg('k')
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', 'k');

f1 = MB_Analysis(pathTofolder,f1);
f1 = NeuronsPhase_Analysis(pathTofolder,f1);
f1 = DNS_Parameters_analysis(pathTofolder,f1);

function f1 = MB_Analysis(pathTofolder,f1)

load(fullfile(pathTofolder,'DataMovingBump'));


Ttask = 4;
plottin = 1;
TeoDur = [450 850 450 850];
TSerialOrder = 4;
TTActivations = TActivations;


%%%%Get AP order across task and serial order
for task = 1:Ttask

    clear  Neuidx1 Tpeak
    Neuidx = IndexNeurons4Condi{task};
    TTneurons = numel(Neuidx);

    for Serial_order = 1:TSerialOrder
        for cellidx = 1:TTneurons

            cellidxOK = Neuidx(cellidx);
            numac = TTActivations{cellidxOK,task}{1,Serial_order}.MaxNumAc;
            if numac
                Tpeak(cellidx,Serial_order) = TTActivations{cellidxOK,task}{1,Serial_order}.PeakTime(1);
            else
                Tpeak(cellidx,Serial_order) = NaN;
            end
        end
    end
    [SortedCells{task},CellGroups{task},PeakQuarterNumAll{task}] = getSortedCellsALLSO(Tpeak,TapTimes(task,:),task,Serial_order);

end



for task = 1:Ttask

    clear  temAcOthers  Neuidx1
    Neuidx = IndexNeurons4Condi{task};
    TTneurons = numel(Neuidx);

    ColorCells = 0;

    taskaligment = 0;%getTaskAligment(DurModCellGroup,task);
    AllPeakQuarterNum = [];

    for Serial_order = 1:TSerialOrder

        [TimeSO,sdinit,sdend] = getsdfTimesSO(task,Serial_order,BehTimesSDF);

        clear Tbegin TDUr TDischargeRate Tpeak PeakMag CellIds SDFCell WinnerSensory AllPeaksSDFtrial
        for cellidx = 1:TTneurons

            cellidxOK = Neuidx(cellidx);
            numac = TTActivations{cellidxOK,task}{1,Serial_order}.MaxNumAc;
            %   cellidx
            if numac
                Tbegin(cellidx,1) = TTActivations{cellidxOK,task}{1,Serial_order}.ActBeg(numac);
                TDUr(cellidx,1) = TTActivations{cellidxOK,task}{1,Serial_order}.ActDur(numac);
                TDischargeRate(cellidx,1)  = TTActivations{cellidxOK,task}{1,Serial_order}.Magni(numac);
                Tpeak(cellidx,1) = TTActivations{cellidxOK,task}{1,Serial_order}.PeakTime(1);
                PeakMag(cellidx,1) = TTActivations{cellidxOK,task}{1,Serial_order}.PeakMag(1);
                TotalDischargeRate(cellidx,1)  = TTActivations{cellidxOK,task}{1,Serial_order}.OverallDischargeR(1);
                AllPeaksSDFtrial(cellidx,1:4)  = TTActivations{cellidxOK,task}{1,Serial_order}.AllPeaksSDFtrial;
                ALLFanoFactor(cellidx,1:4)  = TTActivations{cellidxOK,task}{1,Serial_order}.AllFanoFactor;
                TF = isfield(TTActivations{cellidxOK,task}{1,Serial_order},'TWinnerSensory');
                if TF %sensory poisson
                    WinnerSensory(cellidx,1)  = TTActivations{cellidxOK,task}{1,Serial_order}.TWinnerSensory ;
                else
                    WinnerSensory(cellidx,1)  = 0;
                end


            else
                Tbegin(cellidx,1) = 99999;
                TDUr(cellidx,1) = 99999;
                TDischargeRate(cellidx,1)  = 99999;
                Tpeak(cellidx,1) = 99999;
                PeakMag(cellidx,1) = 99999;
                TotalDischargeRate(cellidx,1)  = 99999;
                AllPeaksSDFtrial(cellidx,1:4)  = 99999;
                ALLFanoFactor(cellidx,1:4)  = 99999;
                WinnerSensory(cellidx,1)  = 99999;
            end
        end

        %   taskaligment = 2;
        if taskaligment == 0 %no aligment
            [out,idxSO]  = sort(Tpeak);
            ColorCellGroups = ones(sum(out < 99990),1);
            temPeak = Tpeak;
            temPeak(temPeak==99999) = [];
            tCellIds = 1:1:numel(temPeak);
            PeakQuarterNum = getPeakQuarterNum3(tCellIds,temPeak,task,Serial_order,TapTimes(task,:));

            Tquarter = 4;
        else
            mColorCellGroups = CellGroups{taskaligment};
            idxSO = SortedCells{taskaligment};
            out = idxSO;
            ColorCellGroups = mColorCellGroups(idxSO);
            PeakQuarterNum =  PeakQuarterNumAll{taskaligment};
            Tquarter = 6;
        end

        CellIds(:,1) = 1:1:TTneurons;

        CellInxfrom2045 = Neuidx(idxSO);
        Sbegin =Tbegin(idxSO);
        SDUr =TDUr(idxSO);
        SDischargeRate =TDischargeRate(idxSO);
        Speak =Tpeak(idxSO);
        SPeakMag = PeakMag(idxSO);
        SCellIds = CellIds(idxSO);
        STotalDischargeRate = TotalDischargeRate(idxSO);
        SAllPeaksSDFtrial = AllPeaksSDFtrial(idxSO,:);
        SALLFanoFactor = ALLFanoFactor(idxSO,:);
        SWinnerSensory = WinnerSensory(idxSO,:);
        AllNeu = 1:numel(SALLFanoFactor(:,1));


        CellInxfrom2045(out == 99999) = [];
        Sbegin(out == 99999) = [];
        SDUr(out == 99999) = [];
        SDischargeRate(out == 99999) = [];
        Speak(out == 99999) = [];
        SPeakMag(out == 99999) = [];
        SCellIds(out == 99999) = [];
        STotalDischargeRate(out == 99999) = [];

        SAllPeaksSDFtrial(out == 99999,:) = [];
        SALLFanoFactor(out == 99999,:) = [];
        SWinnerSensory(out == 99999,:) = [];
        out(out == 99999) = [];

        [Inter,a,b] = intersect(CellInxfrom2045,ColorCells,'rows');
        CNeuidx  = a;


        temPeak = Speak;
        temPeak(temPeak==99999) = [];
        tCellIds = 1:1:numel(temPeak);
        PeakQuarterNumPie = getPeakQuarterNum3(tCellIds,temPeak,task,Serial_order,TapTimes(task,:));
        AllPeakQuarterNum = [AllPeakQuarterNum; PeakQuarterNumPie];

        if plottin
            if task < 3 %for auditory condition
                %%%For autidoty plots of the paper:
                f1 = PlotMovingBumpsJul2023Aud(Sbegin,SDUr,SDischargeRate,Speak,SPeakMag,task,Serial_order,TapTimes(task,1:5),PeakQuarterNum,f1);
            end

        end


        if taskaligment > 0 %no aligment
            out = Speak;
            out(ColorCellGroups>1) =99999;
            CellInxfrom2045(out == 99999) = [];
            Sbegin(out == 99999) = [];
            SDUr(out == 99999) = [];
            SDischargeRate(out == 99999) = [];
            Speak(out == 99999) = [];
            SPeakMag(out == 99999) = [];
            SCellIds(out == 99999) = [];
            STotalDischargeRate(out == 99999) = [];

            SAllPeaksSDFtrial(out == 99999,:) = [];
            SALLFanoFactor(out == 99999,:) = [];
            SWinnerSensory(out == 99999,:) = [];
            SSDFCell(out == 99999,:) = [];
            out(out == 99999) = [];
        end

        MB.TotalActivatedCells(task,Serial_order) = numel(out);
        MB.ActivationDurationM(task,Serial_order) = mean(SDUr);
        MB.ActivationDurationSD(task,Serial_order) = std(SDUr);
        MB.ActivationMagM(task,Serial_order) = mean(SDischargeRate);
        MB.ActivationMagSD(task,Serial_order) = std(SDischargeRate);
        MB.PeakMagM(task,Serial_order) = mean(SPeakMag);
        MB.PeakMagSD(task,Serial_order) = std(SPeakMag);
        MB.TotalDischargeRateM(task,Serial_order) = mean(STotalDischargeRate);
        MB.TotalDischargeRateSD(task,Serial_order) = std(STotalDischargeRate);
        MB.TapTimes(task,Serial_order) = {TapTimes(task,:)}
        MB.CellInxfrom2045(task,Serial_order) = {CellInxfrom2045};


        ploing = 0;
        temRate = ComputeActivationRateMay2021(Speak,ploing);

        ploing = 0;
        PlotOrder = 1;
        [OutCircStatistics] = ComputePeakPolarDistri(Speak,TeoDur(task),temRate.PeakRateBin,temRate.NumThresh75,task,ploing,PlotOrder);
        %[OutCircStatistics] = ComputePeakPolarDistri(Speak,TeoDur(task),temRate.PeakRateBin,temRate.NumThresh75,task,Serial_order,ploing,PlotOrder,TapTimes(task,Serial_order));
        PlotOrder = 10;
        [OutPhasePeakMag] = ComputePhasePeakParam(Speak,SPeakMag,TeoDur(task),task,Serial_order,ploing,PlotOrder,'Peak Magnitude',TapTimes(task,Serial_order));
        PlotOrder = 15;
        [OutPhasePeakDur] = ComputePhasePeakParam(Speak,SDUr,TeoDur(task),task,Serial_order,ploing,PlotOrder,'Peak Duration',TapTimes(task,Serial_order));

        PeakSkew = ComputePeakSkewAPDur(Speak,Sbegin,SDUr);
        PlotOrder = 20;
        [OutPhasePeakSkew] = ComputePhasePeakParam(Speak,PeakSkew,TeoDur(task),task,Serial_order,ploing,PlotOrder,'Peak Skeewness',TapTimes(task,Serial_order));

        PlotOrder = 25;
        [OutPhasePeakTimeVar] = ComputePhasePeakParam(Speak,SAllPeaksSDFtrial(:,2),TeoDur(task),task,Serial_order,ploing,PlotOrder,'Peak Time STD',TapTimes(task,Serial_order));
        PlotOrder = 30;
        [OutPhasePeakMagVar] = ComputePhasePeakParam(Speak,SAllPeaksSDFtrial(:,4),TeoDur(task),task,Serial_order,ploing,PlotOrder,'Peak Mag STD',TapTimes(task,Serial_order));

        PlotOrder = 35;
        [OutPhaseFanoFactor] = ComputePhasePeakParam(Speak,SALLFanoFactor(:,1),TeoDur(task),task,Serial_order,ploing,PlotOrder,'Fano Factor AP',TapTimes(task,Serial_order));

        PlotOrder = 40;
        OutCircStatisticsPolarAP = ComputePolarActivationPeriod(Sbegin,SDUr,TeoDur(task),task,Serial_order,TapTimes(task,Serial_order),ploing,PlotOrder);

        PlotOrder = 50;
        [OutCtePhase] = ComputeCtePhaseParam(Speak,SPeakMag,TeoDur(task),task,Serial_order,ploing,PlotOrder,'Phase CteNumber');

        PlotOrder = 55;
        [OutPhasePeakRateRecritment] = ComputePhasePeakParam(Speak(1:end-1),temRate.Rate,TeoDur(task),task,Serial_order,ploing,PlotOrder,'Cell Recuitment Rate',TapTimes(task,Serial_order));


        MB.MeanRateActivation(task,Serial_order) = temRate.MRate;
        MB.PeakBinRateActivation(task,Serial_order) = temRate.PeakRateBin;
        MB.PeakRateMag(task,Serial_order) = temRate.PeakRateMag;
        MB.WeithedRateMean(task,Serial_order) = temRate.WeithedRateMean;
        MB.MiddleX(task,Serial_order) = temRate.MiddleX;
        MB.LG4PSE(task,Serial_order) = temRate.LG4PSE;
        MB.LG4MaxAsyn(task,Serial_order) = temRate.LG4MaxAsyn;

        MB.NumThresh25(task,Serial_order)  = temRate.NumThresh25;
        MB.NumThresh75(task,Serial_order)  = temRate.NumThresh75;
        MB.NumThresh50(task,Serial_order)  = temRate.NumThresh50;

        MB.TotalCellsQuadrants(task,Serial_order)  = {temRate.totalCells};

        MB.RateInfo{task,Serial_order} = temRate;
        MB.Cellids{task,Serial_order} = SCellIds;
        MB.PeakTimes{task,Serial_order} = Speak;
        MB.PredictedPeakTimes{task,Serial_order} = temRate.ypred{1};

        MB.PeakPhaseMeanResultant(task,Serial_order)  = {OutCircStatistics};
        MB.PolarAPMeanResultant(task,Serial_order)  = {OutCircStatisticsPolarAP};
        MB.PeakPhase(task,Serial_order)  = {OutPhasePeakMag};
        MB.OutPhasePeakDur(task,Serial_order)  = {OutPhasePeakDur};
        MB.OutPhasePeakSkew(task,Serial_order)  = {OutPhasePeakSkew};
        MB.OutPhasePeakTimeVar(task,Serial_order)  = {OutPhasePeakTimeVar};
        MB.OutPhasePeakMagVar(task,Serial_order)  = {OutPhasePeakMagVar};
        MB.OutPhaseFanoFactor(task,Serial_order)  = {OutPhaseFanoFactor};
        MB.OutCteStepPhase(task,Serial_order)  = {OutCtePhase};
        MB.OutPhasePeakRateRecritment(task,Serial_order)  = {OutPhasePeakRateRecritment};

        MB.All.ActivationDurationM{task,Serial_order} = (SDUr);
        MB.All.ActivationBeginM{task,Serial_order} =Sbegin;
        MB.All.ActivationMagM{task,Serial_order} = (SDischargeRate);
        MB.All.TotalDischargeRateM{task,Serial_order} = STotalDischargeRate;
        MB.WinnerSensory{task,Serial_order} = SWinnerSensory;


        clear temRate temPolarActivation
        clear OutCircStatistics OutPhasePeakMag OutPhasePeakDur PeakSkew OutPhasePeakSkew OutPhasePeakTimeVar
        clear OutPhasePeakMagVar OutPhaseFanoFactor  OutCircStatisticsPolarAP

    end

    MB.AllPeakQuarterNum = AllPeakQuarterNum;

    % figure
    plotPieQuarterNum(AllPeakQuarterNum);
end

MBAllSig12 = MB;
return;


function f1 = NeuronsPhase_Analysis(pathTofolder,f1)
%%PlotNeuronsPhaseALLSO2023.m
%%Panel B Figure 5
%%ploting the number of active cells for each phase in a polar plot
%%Hugo Merchant Sept 5/2023
load(fullfile(pathTofolder,'MovinBumpSig'));


idxPlot = [3,4,7,8];

plot_rows = 5;%1
plot_cols = 4;%2

Ttask = 4;
TSerialOrder = 4;
tInterval = [450 850 450 850];
MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';


% clf(figure(4))
colormap('summer');
cmap = colormap;
ncmap = numel(cmap(:,1));
% whitebg([0 0 0])
ColorsAngle = [1:round(ncmap)/16:ncmap];
quartercolors = [1 40 120 240 180 1];



MBAllSig = MBAllSig12;


Modality = [1 1 2 2];
Duration = [450 850 450 850];
m=1;
l = 1;
set(0,'CurrentFigure',f1)


for task = 1:Ttask
    timestep = tInterval(task)/16;

    n = 1;
    t = timestep;
    seqPI = 0;
    clear Angles PhaseM PeakPhaseM PeakPhaseSD binAcRate PeakPhaseMC PhaseSO  PeakPhaseMag NumCellsperAngles BinnedTime
    clear PeakPhaseDur PeakSkew PeakTimeSD_sdf PeakMagSD_sdf FanoFactorPeakPhase RateRecritment PeakPhaseMSO
    for Serial_order = 1:TSerialOrder

        Angles(Serial_order,:) = MBAllSig.PeakPhase{task,Serial_order}(:,2);
        NumCellsperAngles(Serial_order,:) = MBAllSig.PeakPhase{task,Serial_order}(:,7);
        RateRecritment(Serial_order,:) = MBAllSig.OutPhasePeakRateRecritment{task,Serial_order}(:,8);


        t = (tInterval(task)*Serial_order)+timestep;
        seqPI = seqPI+2*pi;
    end


    xp = [0 (2*pi) (4*pi) (6*pi) (8*pi);0 (2*pi) (4*pi) (6*pi) (8*pi) ];
    yp = [0 0 0 0 0; 100 100 100 100 100];

    NumCellsperAngles(4,16) = NumCellsperAngles(3,16);
    tAngles = mean(Angles);
    tNumCellsperAngles = mean(NumCellsperAngles);
    tRateRecritment = mean(RateRecritment);


    subplot(plot_rows,plot_cols,idxPlot(task));

    for ii = 1:16

        r_anglet = tAngles(ii);
        resultantt = tNumCellsperAngles(ii);

        theta = [r_anglet r_anglet];
        rho = [0 resultantt];

        polarplot(theta,rho,'LineWidth', 2, 'Color', [cmap(ColorsAngle(ii),:)]);hold on

    end
    set(gca,'color','k')
    ax(task) = gca;
    ax(task).ThetaColor = 'w';
    ax(task).RColor = 'w';
    ax(task).GridColor = 'w';
    ax(task).GridAlpha = 0.5;

    maxnum = 60;%110 round(max(tNumCellsperAngles));
    tiscnum = [0 maxnum/3 2*maxnum/3];
    tiscnum =round(tiscnum);
    labelstiscnum = num2str(tiscnum);
    haldmax = maxnum/2;
    title([MIL(task,:)],'FontSize',12);
    %title([MIL(task,:)]);
    pax = gca;
    pax.FontSize = 9;
    pax.LineWidth = 1;
    angles = 0:45:360;
    pax.ThetaTick = angles;
    labels = {'0','','90','','180', '', '-90', ''};
    pax.ThetaTickLabel = labels;
    rlim([0 maxnum]);  %rlim([0 150])
    rticks(tiscnum);
    rticklabels({tiscnum(1) tiscnum(2) tiscnum(3) })
    hold off


end

function f1 = DNS_Parameters_analysis(pathTofolder,f1)
%%Panels C-H Figure 5
%%plotting of the pameters of the activations periods of neural sequences
%%of the four task conditions as a function of interval phase
%%Computing ANOVAS for each parameter on duration, modality and interval
%%quarter

%%Hugo Merchant Sept 5/2023
load(fullfile(pathTofolder,'MovinBumpSig'));


% idxPlot = [9,10];
plot_rows = 5;%1
plot_cols = 4;%2

Ttask = 4;
TSerialOrder = 4;
tInterval = [450 850 450 850];
MIL(1,:) = 'A450';
MIL(2,:) = 'A850';
MIL(3,:) = 'V450';
MIL(4,:) = 'V850';

% colormap('jet');
% cmap = colormap;
cmap = [0,0.7490,1; 0 0 1; 1 .5020 0;1 0 0];
ColorTask =[100 10 190 250];
%ColorTask =[20 10 45 60];
% close all;
% whitebg([1 1 1])

xp = [0 (2*pi) (4*pi) (6*pi) (8*pi);0 (2*pi) (4*pi) (6*pi) (8*pi) ];
xp2 = [5 21;5 21];
tModality = [1 1 2 2];
tDuration = [1 2 1 2];

MBAllSig = MBAllSig12;

[AnovaMatrix,RelTimeMatrix] = ComputeAnovaMatrix(MBAllSig);

tTapTimes = MBAllSig.TapTimes;

clear Inter_SDTimePeak InterFanoFactor Inter_APMag

clear MedianPeakPhase
TotalPhaseVal = 64;
TotalInterpolatedPoints = TotalPhaseVal*10;
tbinSO = TotalInterpolatedPoints/TSerialOrder;
tbinSO16 = TotalInterpolatedPoints/(TSerialOrder*10);

m=1;
l = 1;
u = 1;
for task = 1:Ttask
    timestep = tInterval(task)/16;
    n = 1;
    t = timestep;
    xt = [0 (1*tInterval(2)) (2*tInterval(2)) (3*tInterval(2)) (4*tInterval(2));...
        0 (1*tInterval(2)) (2*tInterval(2)) (3*tInterval(2)) (4*tInterval(2))];
    
    initphase = mean(diff(AnovaMatrix(l:l+63,5)));
    PhaseM = AnovaMatrix(l:l+63,5)+initphase;
    PeakTimeSD_sdf = AnovaMatrix(l:l+63,7);%Peak SDdeviation in time
    PeakMagSD_sdf = AnovaMatrix(l:l+63,8);%Peak SDdeviation in magnitude
    FanoFactorPeakPhase = AnovaMatrix(l:l+63,9);%Fano Factor
    PeakPhaseMag = AnovaMatrix(l:l+63,10);%AP magnitude
    NumCellsperAngles = AnovaMatrix(l:l+63,11);%number of cells
    PeakPhaseDur = AnovaMatrix(l:l+63,12); %AP DUration
    PeakSkew = AnovaMatrix(l:l+63,13);%Peak skew withing AP
    RateRecritment = AnovaMatrix(l:l+63,14);%Cell Recruitment rate
    Time = AnovaMatrix(l:l+63,17);%Time Bin
    l = l+TotalPhaseVal;
    dospi = PhaseM(16);
    
    yp = [0 0 0 0 0 ; 120 120 120 120 120];
    MedrateAll = medfilt1(PeakTimeSD_sdf,1); %median filter order 2
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    PeakTimeSD_sdf2 = SmoothRate;
    
    % figure(1)
    set(0,'CurrentFigure',f1)
    ax(9) = subplot(plot_rows,plot_cols,11);
    PeakTimeSD_sdf2(PeakTimeSD_sdf2>95) = 95;
    hold on
    plot(xp,yp,'-w');
    plot(PhaseM,PeakTimeSD_sdf2,'-','LineWidth', 1.5, 'Color', [cmap(task,:)]);

    
    set(ax(9), 'XTick', [xp(1,:)])
    set(ax(9), 'YTick', [40:20:100])
    set(ax(9), 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    set(ax(9), 'YTickLabel', [40:20:100])
    axis([0 26 35 95]);
    % axis square;
    % title(['STD Peak Time' MIL(task,:)]);
    set(ax(9),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(9),'FontSize',13,'LineWidth', 1.5)
    ylabel('SD Peak Time (ms)')
    xlabel('Phase')
    hold off
    % 
    %%%%%Activation rates Interpolated at 440 cells for all durations and SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;    
    for Serial_order = 1:TSerialOrder        
        tem(1:tbinSO16,Serial_order) = PeakTimeSD_sdf2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = PeakTimeSD_sdf2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_SDTimePeak{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) = PeakTimeSD_sdf2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) = PeakTimeSD_sdf2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) = PeakTimeSD_sdf2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) = PeakTimeSD_sdf2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) = PeakTimeSD_sdf2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    
    
    Inter_SDTimePeak2{task} = tem;
    Inter_SDTimePeakM = nanmean(Inter_SDTimePeak2{task}');
    Inter_SDTimePeakSD = nanstd(Inter_SDTimePeak2{task}')/sqrt(16);
    
    % figure(11)
    set(0,'CurrentFigure',f1)
    ax(10) = subplot(plot_rows,plot_cols,12);
    hold on
    Inter_PeakM = mean(Inter_SDTimePeak{task}');
    yp2 = [0 0 ; 120 120];
    xx = 1:1:tbinSO16+8;
    boundedlineHM(ax(10),xx,Inter_SDTimePeakM,Inter_SDTimePeakSD,'cmap',[cmap(task,:)],'transparency',.8)
    % errorbar(xx,Inter_SDTimePeakM,Inter_SDTimePeakSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
%     plot(xx,Inter_SDTimePeakM,'-','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5)
%     plot(xp2,yp2,'-w','LineWidth', 1.5);
    % ylabel('SD Peak Time (ms)')
    xlabel('Phase')
    set(ax(10), 'XTick', [5 13 21])
    set(ax(10), 'YTick', [40:20:100])
    set(ax(10), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    set(ax(10), 'YTickLabel', [])
    set(ax(10),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(10),'FontSize',13,'LineWidth', 1.5)
    axis([0 25 35 95]);
    
    hold off
    
    
    MedrateAll = medfilt1(PeakMagSD_sdf,1); %median filter order 2
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    PeakMagSD_sdf2 = SmoothRate;
    % 
    % figure(2)
    % whitebg([0 0 0])    
    % hold on
    % plot(xp,yp,'-w');
    % plot(PhaseM,PeakMagSD_sdf2,'-','LineWidth', 1.5, 'Color', [cmap(ColorTask(task),:)]);
    % 
    % set(gca, 'XTick', [xp(1,:)])
    % set(gca, 'YTick', [0:5:15])
    % set(gca, 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    % set(gca, 'YTickLabel', [0:5:15])
    % axis([0 26 5 15]);    
    % % title(['STD Peak Time' MIL(task,:)]);
    % set(gca,'TickDir','out','TickLength', [0.02 0.02])
    % set(gca,'FontSize',13,'LineWidth', 1.5)
    % ylabel('SD Peak Magnitude (Hz)')
    % xlabel('Phase')
    % hold off
    
    %%%%%Activation rates Interpolated at 440 cells for all durations and SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
    % teminter = cellfun(@InterpRateActivationT, {PeakMagSD_sdf},{TotalInterpolatedPoints} , 'UniformOutput', 0);
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = PeakMagSD_sdf2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = PeakMagSD_sdf2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_SDMagPeak{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) = PeakMagSD_sdf2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) = PeakMagSD_sdf2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) = PeakMagSD_sdf2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) = PeakMagSD_sdf2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) = PeakMagSD_sdf2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    
    Inter_SDMagPeakM2{task} = tem;
    
    %Inter_SDMagPeak{task} = downsample(tem,10);
    Inter_SDMagPeakM = nanmean(Inter_SDMagPeakM2{task}');
    Inter_SDMagPeakSD = nanstd(Inter_SDMagPeakM2{task}')/sqrt(16);
    
    
    % figure(12)
    % whitebg([0 0 0])    
    % hold on
    % %Inter_PeakM = mean(Inter_SDMagPeak{task}');
    % 
    % xx = 1:1:tbinSO16+8;
    % boundedlineHM(xx,Inter_SDMagPeakM,Inter_SDMagPeakSD,'cmap',[cmap(ColorTask(task),:)],'transparency',0.6)
    % % errorbar(xx,Inter_SDTimePeakM,Inter_SDTimePeakSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
    % plot(xx,Inter_SDMagPeakM,'-','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5)
    % ylabel('SD Peak Magnitude (Hz)')
    % xlabel('Phase')
    % plot(xp2,yp2,'-w','LineWidth', 1.5);
    % set(gca, 'XTick', [5 13 21])
    % set(gca, 'YTick', [0:5:15])
    % set(gca, 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    % set(gca, 'YTickLabel', [0:5:15])
    % set(gca,'TickDir','out','TickLength', [0.02 0.02])
    % set(gca,'FontSize',13,'LineWidth', 1.5)
    % axis([0 25 5 15]);
    % 
    % hold off
    % 
    
    MedrateAll = medfilt1(FanoFactorPeakPhase,1); %median filter order 1 whichmean there is NO filter
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    FanoFactorPeakPhase2 = SmoothRate;
    
    
    % figure(3)    
    % hold on
    % yp = [-2 -2 -2 -2 -2 ; 2 2 2 2 2 ];
    % hold on
    % plot(xp,yp,'-w');
    % plot(PhaseM,FanoFactorPeakPhase2,'-','LineWidth', 1.5, 'Color', [cmap(ColorTask(task),:)]);
    % 
    % set(gca, 'XTick', [xp(1,:)])
    % set(gca, 'YTick', [0:0.5:2])
    % set(gca, 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    % set(gca, 'YTickLabel', [0:0.5:2])
    % axis([0 26 1 2]);
    % % title(['Fano Factor ' MIL(task,:)]);
    % set(gca,'TickDir','out','TickLength', [0.02 0.02])
    % set(gca,'FontSize',13,'LineWidth', 1.5)
    % ylabel('Fano factor')
    % xlabel('Phase')
    % hold off
    % 
    %%%%%FanoFactorPeakPhase mean four SO,
    %%%%%plotted for task
    %%%%%Activation rates Interpolated at 440 cells for all durations and SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
    % teminter = cellfun(@InterpRateActivationT, {PeakMagSD_sdf},{TotalInterpolatedPoints} , 'UniformOutput', 0);
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = FanoFactorPeakPhase2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = FanoFactorPeakPhase2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_FanoFactor{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) = FanoFactorPeakPhase2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) = FanoFactorPeakPhase2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) = FanoFactorPeakPhase2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) = FanoFactorPeakPhase2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) = FanoFactorPeakPhase2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    
    Inter_FanoFactor2{task} = tem;
    Inter_FanoFactorM = nanmean(Inter_FanoFactor2{task}');
    Inter_FanoFactorSD = nanstd(Inter_FanoFactor2{task}')/sqrt(16);
    
    
    %figure(13)
    set(0,'CurrentFigure',f1)
    ax(12) = subplot(plot_rows,plot_cols,15);
    hold on
    Inter_PeakM = mean(Inter_Peak');
    
    xx = 1:1:tbinSO16+8;
    boundedlineHM(ax(12),xx,Inter_FanoFactorM,Inter_FanoFactorSD,'cmap',[cmap(task,:)],'transparency',0.6)
    %errorbar(xx,Inter_FanoFactorM,Inter_FanoFactorSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
    plot(xx,Inter_FanoFactorM,'-','Color',[cmap(task,:)],'LineWidth',1.5)
    %  plot(xx, Inter_PeakM,'-k');
    
    ylabel('Fano factor')
    xlabel('Phase')
    plot(xp2,yp2,'-w','LineWidth', 1.5);
    set(ax(12), 'XTick', [5 13 21])
    set(ax(12), 'YTick', [1.2:0.2:1.8])
    set(ax(12), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    set(ax(12), 'YTickLabel', [1.2:0.2:1.8])
    set(ax(12),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(12),'FontSize',13,'LineWidth', 1.5)
    axis([0 25 1.2 1.8]);
    
    hold off
    
    MedrateAll = medfilt1(PeakPhaseMag,1); %median filter order 2
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    PeakPhaseMag2 = SmoothRate;
    
    %%%%%%%%%AP Discharge rat
    % figure(4)
    % %subplot(1,2,1)
    % hold on
    % yp = [0 0 0 0 0 ; 35 35 35 35 35];
    % PeakPhaseMag2(PeakPhaseMag2>32) = 32;
    % plot(xp,yp,'-w','LineWidth', 1.5);
    % plot(PhaseM,PeakPhaseMag2,'-','LineWidth', 1.5, 'Color', [cmap(ColorTask(task),:)]);
    % set(gca, 'XTick', [xp(1,:)])
    % set(gca, 'YTick', [0:5:35])
    % set(gca, 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    % set(gca, 'YTickLabel', [0:5:35])
    % axis([0 26 5 32]);
    % % title(['Discharge Rate Per Angle ' MIL(task,:)]);
    % set(gca,'TickDir','out','TickLength', [0.02 0.02])
    % set(gca,'FontSize',13,'LineWidth', 1.5)
    % ylabel('AP Discharge rate')
    % xlabel('Phase')
    % hold off
    % xpp = [0 2.5];
    % ypp = [0 0];
    
    %%%%%Magnitude Activation period mean four SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
    % teminter = cellfun(@InterpRateActivationT, {PeakMagSD_sdf},{TotalInterpolatedPoints} , 'UniformOutput', 0);
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = PeakPhaseMag2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = PeakPhaseMag2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_APMag{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) = PeakPhaseMag2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) = PeakPhaseMag2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) = PeakPhaseMag2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) = PeakPhaseMag2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) = PeakPhaseMag2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    Inter_APMag2{task} = tem;
    Inter_APMagM = nanmean(Inter_APMag2{task}');
    Inter_APMagSD = nanstd(Inter_APMag2{task}')/sqrt(16);
    
    
    % figure(14) 
    set(0,'CurrentFigure',f1)
    ax(11) = subplot(plot_rows,plot_cols,14);
    %
    hold on
    tInter_APMagM = mean(Inter_Peak');
    
    xx = 1:1:tbinSO16+8;
    boundedlineHM(ax(11),xx,Inter_APMagM,Inter_APMagSD,'cmap',[cmap(task,:)],'transparency',0.8)
    %errorbar(xx,Inter_APMagM,Inter_APMagSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1);
    plot(xx,Inter_APMagM,'-','Color',[cmap(task,:)],'LineWidth',1)
    %  plot(xx, Inter_PeakM,'-k');
    
    ylabel('AP Discharge rate')
    xlabel('Phase')
    plot(xp2,yp2,'-w','LineWidth', 1.5);
    set(ax(11), 'XTick', [5 13 21])
    set(ax(11), 'YTick', [0:5:35])
    set(ax(11), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    set(ax(11), 'YTickLabel', [0:5:35])
    set(ax(11),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(11),'FontSize',13,'LineWidth', 1.5)
    axis([0 25 5 32]);
    hold off
    
    
    MedrateAll = medfilt1(PeakPhaseDur,1); %median filter order 2
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    PeakPhaseDur2 = SmoothRate;
    
    
    %%%%%Activatio period DURATION
    % figure(5)
    set(0,'CurrentFigure',f1)
    ax(7) = subplot(plot_rows,plot_cols,9);
    hold on
    yp = [100 100 100 100 100 ; 500 500 500 500 500];
    
    plot(ax(7),PhaseM,PeakPhaseDur2,'-','LineWidth', 1.5, 'Color', [cmap(task,:)]);
    plot(ax(7),xp,yp,'-w');    
    set(ax(7), 'XTick', [xp(1,:)])
    set(ax(7), 'YTick', [150:100:550])
    set(ax(7), 'XTickLabel', [{'0'; '2\pi'; '4\pi'; '6\pi'; '8\pi'}])
    set(ax(7), 'YTickLabel', [150:100:550])
    
    axis([0 26 150 450]);    
    % title(['AP Duration Per Angle ' MIL(task,:)]);
    set(gca,'TickDir','out','TickLength', [0.02 0.02])
    set(gca,'FontSize',13,'LineWidth', 1.5)
    ylabel('AP Duration (ms)')
    xlabel('Phase')
    hold off
    xpp = [0 2.5];
    ypp = [0 0];
    
    %%%%%Activatio period DURATION period mean four SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
    % teminter = cellfun(@InterpRateActivationT, {PeakMagSD_sdf},{TotalInterpolatedPoints} , 'UniformOutput', 0);
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = PeakPhaseDur2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = PeakPhaseDur2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_APDur{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) =  PeakPhaseDur2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) =  PeakPhaseDur2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) =  PeakPhaseDur2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) =  PeakPhaseDur2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) =  PeakPhaseDur2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    Inter_APDur2{task} = tem;
    Inter_APDurM = nanmean(Inter_APDur2{task}');
    Inter_APDurSD = nanstd(Inter_APDur2{task}')/sqrt(16);
    
    
    % figure(15)
    set(0,'CurrentFigure',f1)
    ax(8) = subplot(plot_rows,plot_cols,10);
    yp2 = [0 0 ; 500 500];    
    hold on
    tInter_APDurM = mean(Inter_Peak');
    
    xx = 1:1:tbinSO16+8;
    boundedlineHM(ax(11),xx,Inter_APDurM,Inter_APDurSD,'cmap',[cmap(task,:)],'transparency',0.8)
    %errorbar(xx,Inter_APDurM,Inter_APDurSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
    plot(xx,Inter_APDurM,'-','Color',[cmap(task,:)],'LineWidth',1.5)
    %  plot(xx, Inter_PeakM,'-k');
    
    % ylabel('AP Duration (ms)')
    xlabel('Phase')
    plot(xp2,yp2,'-w','LineWidth', 1.5);
    set(ax(8), 'XTick', [5 13 21])
    set(ax(8), 'YTick', [150:100:550])
    set(ax(8), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    set(ax(8), 'YTickLabel', [])
    set(ax(8),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(8),'FontSize',13,'LineWidth', 1.5)
    axis([0 25 150 500]);    
    hold off
    
    MedrateAll = medfilt1(RateRecritment,1); %median filter order 1
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    RateRecritment2 = SmoothRate;
    
    %%%%%%%%RECRUITMEN RATE
    % figure(6)    
    set(0,'CurrentFigure',f1)
    ax(13) = subplot(plot_rows,plot_cols,17);
    hold on
    yp = [0 0 0 0 0 ; 3 3 3 3 3];
    
    RateRecritment2 =  round(RateRecritment2,5);
    plot(PhaseM,RateRecritment2,'-','LineWidth', 1.5, 'Color', [cmap(task,:)]);
    
    plot(xp,yp,'-w');
    plot(xpp,ypp,'-w');
    set(ax(13), 'XTick', [xp(1,:)])
    set(ax(13), 'YTick', [0:1:4])
    set(ax(13), 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    set(ax(13), 'YTickLabel', [0:1:4])
    
    axis([0 26 0.2 4.5]);
    % title(['Cell Recruitment Rate ' MIL(task,:)]);
    set(ax(13),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(13),'FontSize',13,'LineWidth', 1.5)
    ylabel('Recruitment rate')
    xlabel('Phase')
    hold off
    xpp = [0 2.5];
    ypp = [0 0];
    
    %%%%%Activatio period RECRUITMEN RATE period mean four SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
   
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = RateRecritment2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = RateRecritment2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_RecruiRate{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) =   RateRecritment2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) =   RateRecritment2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) =   RateRecritment2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) =   RateRecritment2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) =   RateRecritment2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    Inter_RecruiRate2{task} = tem;
    Inter_RecruiRateM = nanmean(Inter_RecruiRate2{task}');
    Inter_RecruiRateSD = nanstd(Inter_RecruiRate2{task}')/sqrt(16);
    
    
    % figure(16)
    set(0,'CurrentFigure',f1)
    ax(14) = subplot(plot_rows,plot_cols,18);
 
    hold on
    tInter_RecruiRateM = mean(Inter_Peak');
    
    xx = 1:1:tbinSO16+8;
    boundedlineHM(ax(14),xx,Inter_RecruiRateM,Inter_RecruiRateSD,'cmap',[cmap(task,:)],'transparency',0.8)
    %errorbar(xx,Inter_RecruiRateM,Inter_RecruiRateSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
    plot(xx,Inter_RecruiRateM,'-','Color',[cmap(task,:)],'LineWidth',1.5)
    %  plot(xx, Inter_PeakM,'-k');
    
    % ylabel('Recruitment rate')
    xlabel('Phase')
    plot(xp2,yp2,'-w','LineWidth', 1.5);
    set(ax(14), 'XTick', [5 13 21])    
    set(ax(14), 'YTick', [0:1:6])
    set(ax(14), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])
    set(ax(14), 'YTickLabel', [0:1:6])    
    set(ax(14),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(14),'FontSize',13,'LineWidth', 1.5)  
    axis([0 25 0.2 6.5]);    
    hold off
    
    
    MedrateAll = medfilt1(NumCellsperAngles,1); %median filter order 1 which mean no filter
    
    SmoothRate = smooth(MedrateAll, 0.08,'rloess');
    SmoothRate(SmoothRate<0) = 0;
    if sum(SmoothRate) < 10
        SmoothRate = MedrateAll;
    end
    SmoothRate(SmoothRate<0) = 0;
    NumCellsperAngles2 = SmoothRate;
    
    %%%% Number of cells
    % figure(7)   
    set(0,'CurrentFigure',f1)
    ax(15) = subplot(plot_rows,plot_cols,19);
    hold on
    yp = [0 0 0 0 0; 140 140 140 140 140];
    plot(xp,yp,'-w','LineWidth', 1.5);
    plot(PhaseM,NumCellsperAngles2,'-','LineWidth', 1.5, 'Color', [cmap(task,:)]);
    
    maxNum = max(NumCellsperAngles);
    MaxScale = round(maxNum + maxNum*.1);
    
    set(ax(15) , 'XTick', [xp(1,:)])
    set(ax(15) , 'YTick', [0:10:160])
    set(ax(15) , 'YTickLabel', [0:10:160])
    set(ax(15) , 'XTickLabel', [{'0'; '2*pi'; '4*pi'; '6*pi'; '8*pi'}])
    
    
    axis([0 26 0 56]);
    
    % title(['Number of Cells Per Angle ' MIL(task,:)]);
    set(ax(15) ,'TickDir','out','TickLength', [0.02 0.02])
    set(ax(15) ,'FontSize',13,'LineWidth', 1.5)
    ylabel('Number of cells')
    xlabel('Phase')
    hold off
    
    %%%%%Activatio period number of cells period mean four SO,
    %%%%%plotted for task
    clear  tem
    nn=1;
    nnn=1;
    % teminter = cellfun(@InterpRateActivationT, {PeakMagSD_sdf},{TotalInterpolatedPoints} , 'UniformOutput', 0);
    for Serial_order = 1:TSerialOrder
        tem(1:tbinSO16,Serial_order) = NumCellsperAngles2(nnn:nnn+tbinSO16-1,1);
        Inter_Peak(1:tbinSO16,Serial_order) = NumCellsperAngles2(nnn:nnn+tbinSO16-1,1);
        nn= nn+tbinSO;
        nnn= nnn+tbinSO16;
    end
    Inter_NumCells{task} = tem;
    
    clear  tem
    nn=1;
    nnn=1;
    
    tem(1:4,1) = NaN;
    
    for Serial_order = 1:1%TSerialOrder
        tem(5:tbinSO16+8,Serial_order) =   NumCellsperAngles2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 2:3
        tem(1:4,Serial_order) =  NumCellsperAngles2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+8,Serial_order) =   NumCellsperAngles2(nnn:nnn+tbinSO16+3,1);
        nnn= nnn+tbinSO16;
    end
    for Serial_order = 4:4
        tem(1:4,Serial_order) =   NumCellsperAngles2(nnn-4:nnn-1,1);
        tem(5:tbinSO16+4,Serial_order) =   NumCellsperAngles2(nnn:nnn+tbinSO16-1,1);
    end
    tem(21:24,4) = NaN;
    
    Inter_NumCells2{task} = tem;
    Inter_NumCellsM = nanmean(Inter_NumCells2{task}');
    Inter_NumCellsSD = nanstd(Inter_NumCells2{task}')/sqrt(16);
    
    delta = 7;
    time = (1:1:18)';
  
    model = 5;%gaussian
    Inter_NumCellsMM = Inter_NumCellsM;
    Inter_NumCellsMM(1:6) = [];
    [temcoe,resid,ypred] = curvefitting6functions2019(time,Inter_NumCellsMM',model,task);
    temradm = temcoe(1);%%mean gaussian
    temradsd = temcoe(2);
    
    unitpifor16 = dospi/16;
    radians20units = unitpifor16*20;
    timeseriesradians20 = unitpifor16:unitpifor16:radians20units;
    timeseriesradians18 = timeseriesradians20(3:end)';
    [temcoer,resid,ypredr] = curvefitting6functions2019(timeseriesradians18,Inter_NumCellsMM',model,task);
    
    
    temcoe(6) = temcoer(1);%mean radians
    temcoe(7) = temcoer(2);%sigma radians.
    
    coef(1:7,task) = temcoe;
    
    
    % figure(17)  
    set(0,'CurrentFigure',f1)
    ax(16) = subplot(plot_rows,plot_cols,20);
    hold on;
%     tInter_NumCellsM = mean(Inter_Peak');
%     Inter_NumCellsSD = std(Inter_Peak');
    
    xx = 1:1:tbinSO16+8;
    time = (7:1:24)';
    boundedlineHM(ax(16),xx,Inter_NumCellsM,Inter_NumCellsSD,'cmap',[cmap(task,:)],'transparency',.8)
    %errorbar(xx,Inter_NumCellsM,Inter_NumCellsSD,'.','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5);
    %  plot(xx,Inter_NumCellsM,'-','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5)    
    %    plot(time,ypred,'-','Color',[cmap(ColorTask(task),:)],'LineWidth',1.5)
    %    %%Predicted Gaussian
    %  plot(xx, Inter_PeakM,'-k');
    
    % ylabel('Number of cells')
    xlabel('Phase')
    plot(xp2,yp2,'-w','LineWidth', 1.5);
    set(ax(16), 'XTick', [5 13 21])    
    set(ax(16), 'YTick', [0:10:160])
    set(ax(16), 'XTickLabel', [{'0'; 'pi'; '2*pi'}])    
    set(ax(16), 'YTickLabel', [])
    set(ax(16),'TickDir','out','TickLength', [0.02 0.02])
    set(ax(16),'FontSize',13,'LineWidth', 1.5)    
    
    axis([0 25 0 56]);
    hold off
    
    
    TapTimes = round(tTapTimes{task}(1:5));
    ProducedIntervals(task,:) = diff(TapTimes);
    for so =1 :4
        allpeak = MBAllSig.PeakTimes{task,so};
        tallpeak = allpeak - TapTimes(so);
        tallpeaknext = allpeak - TapTimes(so+1);
        MedianPeakPhase(task,so) = median(tallpeak);
        MedianPeakPhaseNext(task,so) = median(tallpeaknext);
        ModePeakPhaseNext(task,so) = mode(tallpeaknext);
        
        temdata = (tallpeak.*360)./ProducedIntervals(task,so); %Results in deegrees
        temdata(temdata> 359.999999) = 360;
        MedianPeakPhaseDegrees(task,so) = median(temdata);
        valsRad = deg2rad(temdata); %Results in radians
        
        Nbins = 16;
        stepedges = 2*pi/((Nbins));
        edgesPolar = 0:stepedges:2*pi;
        edgesPolar2 = edgesPolar;
        edgesdegrees = rad2deg(edgesPolar);
        N = histc(valsRad,edgesPolar);
        N(16) = N(16)+N(17);
        edgesPolar(end) = [];
        N(end) = [];
        hiscoun = N';
        
    end
    
    
    
    
    
    
end

%New Parameteres
set(ax([8,10,11,12,16]), 'XTickLabel', [{'0'; '\pi';'2\pi'}])
set(ax([7,9,13,15]), 'XTickLabel', [{'0'; '2\pi'; '4\pi'; '6\pi'; '8\pi'}])

set(ax(7:16),'Linewidth',1.5)
set(ax(7:16),'TickDir','out');
set(ax(7:16),'FontSize',10)


%%%%ANOVAS for each parameter

clear anTvar
nn = 1;

cellquarter = TotalPhaseVal/16;

ratesA = Inter_SDTimePeak;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
clear AnovaStatDurModQuarter
[pvalues_MBQuarter(1:6,1),AnovaStatDurModQuarter{1}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');


clear anTvar
nn = 1;

cellquarter = TotalPhaseVal/16;

ratesA = Inter_SDMagPeak;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,2),AnovaStatDurModQuarter{2}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');



clear anTvar
nn = 1;
ratesA = Inter_FanoFactor;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,3),AnovaStatDurModQuarter{3}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');


clear anTvar
nn = 1;
ratesA = Inter_APMag;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,4),AnovaStatDurModQuarter{4}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');

clear anTvar
nn = 1;
ratesA = Inter_APDur;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,5),AnovaStatDurModQuarter{5}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');


clear anTvar
nn = 1;
ratesA = Inter_RecruiRate;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,6),AnovaStatDurModQuarter{6}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');


clear anTvar
nn = 1;
ratesA = Inter_NumCells;

for task = 1:Ttask
    nq = 1;
    for quarter = 1:4
        arquar = mean(ratesA{task}(nq:nq+cellquarter-1,:)');
        anTvar(nn:nn+3,1) = arquar;
        anTvar(nn:nn+3,2) = tModality(task);
        anTvar(nn:nn+3,3) = tDuration(task);
        anTvar(nn:nn+3,4) = quarter;
        nn = nn+4;
        nq =nq+cellquarter;
    end
end
interactionmatrix =  [1 0 0; 0 1 0 ;0 0 1;...
    1 1 0; ; 1 0 1; 0 1 1];
[pvalues_MBQuarter(1:6,7),AnovaStatDurModQuarter{7}] = anovan(anTvar(:,1) ,{anTvar(:,2) anTvar(:,3) anTvar(:,4)},'model',interactionmatrix,'display','off');








