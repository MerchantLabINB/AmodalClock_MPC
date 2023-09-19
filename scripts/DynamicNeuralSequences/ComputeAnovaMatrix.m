function [AnovaMatrix,RelTimeMatrix] = ComputeAnovaMatrix (MBAllSig)


Ttask = 4;
TSerialOrder = 4;
Modality = [1 1 2 2];
tInterval = [450 850 450 850];
nInterval = [1 2 1 2];
m=1;
l = 1;
u = 1;
RelTimeMatrix = [];
for task = 1:Ttask
    timestep = tInterval(task)/16;
    n = 1;
    o = 1;
    t = timestep;
    seqtime = 0;
    seqtimeo = 0;
    acu  = 0;
    seqPI = 0;
    clear TimeM PhaseM PeakPhaseM PeakPhaseSD binAcRate PeakPhaseMC PhaseSO  PeakPhaseMag NumCellsperAngles BinnedTime
    clear PeakPhaseDur PeakSkew PeakTimeSD_sdf PeakMagSD_sdf FanoFactorPeakPhase RateRecritment PeakPhaseMSO
    clear TimeRel850 NumCellsRel850 Binns
    for Serial_order = 1:TSerialOrder
        
        TimeRel850(o:o + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,4)+seqtimeo;
        NumCellsRel850(o:o + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,7);
        
        TimeSO(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,4);
        TimeM(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,4)+seqtime;
        PhaseSO(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,2);
        PhaseM(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,2)+seqPI;
        PeakPhaseM(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,5);
        PeakPhaseMag(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,8);
        NumCellsperAngles(n:n + 15,1) = MBAllSig.PeakPhase{task,Serial_order}(:,7);
        PeakPhaseDur(n:n + 15,1) = MBAllSig.OutPhasePeakDur{task,Serial_order}(:,8);
        PeakSkew(n:n + 15,1) = MBAllSig.OutPhasePeakSkew{task,Serial_order}(:,8);
        
        RateRecritment(n:n + 15,1) = MBAllSig.OutPhasePeakRateRecritment{task,Serial_order}(:,8);
        
        temMC = MBAllSig.PeakPhase{task,Serial_order}(:,5);
        temMC(temMC<0) = (2*pi)+temMC(temMC<0);
        PeakPhaseMC(n:n + 15,1) = temMC+seqPI;
        PeakPhaseMSO(n:n + 15,1) = temMC;
        
        temCte = MBAllSig.OutCteStepPhase{task,Serial_order}(:,3);
        temCte(temCte<0) = (2*pi)+temCte(temCte<0);
        PeakCtePhase(n:n + 15,1) = temCte+seqPI;
        
        PeakTimeSD_sdf(n:n + 15,1) = MBAllSig.OutPhasePeakTimeVar{task,Serial_order}(:,8);
        PeakMagSD_sdf(n:n + 15,1) = MBAllSig.OutPhasePeakMagVar{task,Serial_order}(:,8);
        
        FanoFactorPeakPhase(n:n + 15,1) = MBAllSig.OutPhaseFanoFactor{task,Serial_order}(:,8);
        
        BinnedTime(n:n + 15,1) = MBAllSig.OutPhaseFanoFactor{task,Serial_order}(:,4);
        
        AnovaMatrix(m:m+15,1) = Modality(task);
        AnovaMatrix(m:m+15,2) = tInterval(task);
        AnovaMatrix(m:m+15,3) = Serial_order;
        
        AnovaMatrix(m:m+15,19) = nInterval(task);
       
        Binns(n:n + 15,1)  = 1:16;
        
        m=m+16;
       
        
        binAc =MBAllSig.PeakBinRateActivation(task,Serial_order);
        NumCellsCondition = numel(MBAllSig.PeakTimes{task,Serial_order});
        tembinAcRate = (binAc.*360)./NumCellsCondition; %Results in deegrees
        binAcRate(Serial_order,1) = deg2rad(tembinAcRate)+seqPI; %Results
        
        Num25 = MBAllSig.NumThresh25(task,Serial_order);
        Num50 = MBAllSig.NumThresh50(task,Serial_order);
        Num75 = MBAllSig.NumThresh75(task,Serial_order);
        temNum25 = (Num25.*360)./NumCellsCondition; %Results in deegrees
        temNum50 = (Num50.*360)./NumCellsCondition; %Results in deegrees
        temNum75 = (Num75.*360)./NumCellsCondition; %Results in deegrees
        Phase25(Serial_order,1) = deg2rad(temNum25)+seqPI;
        Phase50(Serial_order,1) = deg2rad(temNum50)+seqPI;
        Phase75(Serial_order,1) = deg2rad(temNum75)+seqPI;
        
        
               
        n = n + 16;
        t = (tInterval(task)*Serial_order)+timestep;
        seqPI = seqPI+2*pi;
        seqtime = (tInterval(task)*Serial_order);
        
        if task == 1 | task == 3
            sp = (tInterval(task)*Serial_order)+acu;
            TimeRel850(o+16) = sp;
            TimeRel850(o+17) = sp +400-timestep;
            NumCellsRel850(o+16) = 0;
            NumCellsRel850(o+17) = 0;
            o = o+18;
            seqtimeo = sp + 400;
            acu = acu+400;
        else
            o = o+16;
            seqtimeo = seqtime;
        end
        
    end
    
    AnovaMatrix(l:l+63,4) = PhaseSO;
    AnovaMatrix(l:l+63,5) = PhaseM;    
    AnovaMatrix(l:l+63,6) = PeakPhaseMSO;%phase
    AnovaMatrix(l:l+63,7) = PeakTimeSD_sdf;%Peak SDdeviation in time
    AnovaMatrix(l:l+63,8) = PeakMagSD_sdf;%Peak SDdeviation in magnitude
    AnovaMatrix(l:l+63,9) = FanoFactorPeakPhase;%Fano Factor
    AnovaMatrix(l:l+63,10) = PeakPhaseMag;%AP magnitude
    AnovaMatrix(l:l+63,11) = NumCellsperAngles;%number of cells
    AnovaMatrix(l:l+63,12) = PeakPhaseDur;%AP DUration
    AnovaMatrix(l:l+63,13) = PeakSkew;%Peak skew withing AP
    AnovaMatrix(l:l+63,14) = RateRecritment;%Cell Recruitment rate
    AnovaMatrix(l:l+63,15) = task;
    AnovaMatrix(l:l+63,16) = BinnedTime;
    AnovaMatrix(l:l+63,17) = TimeM;%time interval/16 across all SO
    AnovaMatrix(l:l+63,18) = Binns;%biins within serial order
    l = l+64;
    
    
   if task == 1 | task == 3
       RelTimeMatrix(u:u+71,1) = TimeRel850;
       RelTimeMatrix(u:u+71,2) = NumCellsRel850;
       u=u+72;
   else
       RelTimeMatrix(u:u+63,1) = TimeRel850;
       RelTimeMatrix(u:u+63,2) = NumCellsRel850;
       u=u+64;
   end
   
end