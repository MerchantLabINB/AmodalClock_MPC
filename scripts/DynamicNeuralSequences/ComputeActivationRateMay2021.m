function ActivationRate = ComputeActivationRateMay2021(Onset,ploing)

%%Onset of Peak

nnMb = numel(Onset);
if  nnMb > 140
    OutTails = 30;
    RemOnset = 10;
elseif  nnMb > 40
   OutTails = 5; 
   RemOnset = 5;
else
    OutTails = 2; 
    RemOnset = 1;
end


colorH = 2;
rateAll = diff(Onset,1);
TimeDiffAllCells = Onset(end) - Onset(1);

MedrateAll = medfilt1(rateAll,2); %median filter order 2

SmoothRate = smooth(MedrateAll, 0.08,'rloess');
SmoothRate(SmoothRate<0) = 0;
if sum(SmoothRate) < 10
    SmoothRate = MedrateAll;
end
SmoothRate(SmoothRate<0) = 0;

xx = 1:1:numel(SmoothRate);

[PeakRate,BinPeakRate] = max(SmoothRate(OutTails-1:end-OutTails));
BinPeakRate = BinPeakRate+(OutTails-2);
if ploing
    figure
    hold on
    plot(xx,MedrateAll,'r-',xx,SmoothRate,'b-');
    plot(BinPeakRate,PeakRate,'.y');
end

% 
% 
% if BinPeakRate > 150
%     [PeakRate,BinPeakRate] = max(SmoothRate(OutTails:150));
%     BinPeakRate = BinPeakRate+OutTails;
%     if ploing
%         plot(BinPeakRate,PeakRate,'.y','MarkerSize', 15);
%     end
% end
% if BinPeakRate < 43
%     [PeakRate,BinPeakRate] = max(SmoothRate(50:140));
%     BinPeakRate = BinPeakRate+50;
%     if ploing
%         plot(BinPeakRate,PeakRate,'.y','MarkerSize', 15);
%     end
% end

%hold off
% SD = normpdf(xx',SmoothRate, STD);
% [m,s] = normfit(xx,SmoothRate(OutTails:end-OutTails));
% [m1,s1] = normfit(xx,SmoothRate*1000);
% ypred = normpdf(xx,m,s);
% ypred1 = normpdf(xx,m1,s1);
% ypred  = ypred *1000;
% plot(xx,ypred,'r-',xx,SmoothRate,'b-');
% plot(xx,ypred1,'g-');
% plot(xx,SD*1000,'r-');



skeqq = skewness(SmoothRate(OutTails:end-OutTails));


xx = (OutTails:1:nnMb-OutTails-1)';
wm = wmean(xx,SmoothRate(OutTails:end-OutTails)*100);

MedOnsetAll = medfilt1(Onset,2);
MedOnsetAll(1) = [];
SmoothOnset = smooth(MedOnsetAll);

MinS = min(SmoothOnset);
SmoothOnset = SmoothOnset-MinS;
SmoothOnset(end-RemOnset:end) = [];
MaxS = max(SmoothOnset);
NorSmoothOnset = SmoothOnset/MaxS*100;
xx0 = (1:1:numel(SmoothOnset))';

PseL = nnMb/2-(nnMb/2*0.4);
PseH = nnMb/2+(nnMb/2*0.4);
Lowerbounds = [-50 0.7 PseL 80];
Higherbounds = [15 2.9 PseH 130];
[cf,G]=L4P(xx0,NorSmoothOnset);
%[cfB,GB]=L4P(xx0,NorSmoothOnset,Lowerbounds,Higherbounds);
ypred = feval(cf,xx0);
%ypredB = feval(cfB,xx0);
%xpred=L4Pinv(cf,NorSmoothOnset');
xpred=(L4Pinv(cf,[25 50 75]));
y255075pred = feval(cf,xpred);
xpred = round(xpred);
ypredround = round(ypred);
numypred25 = xpred(1);%min(find(ypredround >= 24 & ypredround <= 26));
numypred50 = xpred(2);%min(find(ypredround >= 49 & ypredround <= 51));
numypred75 = xpred(3);%min(find(ypredround >= 74 & ypredround <= 76));
totalCells = [numypred25 numypred50-numypred25 numypred75-numypred50 numel(ypredround)-numypred75];


if ploing
    figure
    hold on;
    plot(xx0,NorSmoothOnset,'b-', xx0,ypred,'-g');% xx0,ypredB,'-r');
    plot(numypred25,ypred(numypred25),'ob', numypred75,ypred(numypred75),'ob',numypred50,ypred(numypred50),'og');
    plot(BinPeakRate,ypred(BinPeakRate),'.k','MarkerSize', 15);
end




% 
% 
% Onset0 = Onset-Onset(1);
% [OutCircularStats,OutphaseResults] = 
% [circSTD,r_angle,pvalRalegigh,MeanResultant,outMeanR] = circplot(Onset0,TeoDur,colorH,BinPeakRate,numypred75,condi,SO,ploing)
% MeanR_degrees = rad2deg(r_angle);
% 
% 
% MPeak0 = MPeak-Onset(1);
% [outMeanRMSDPeaksdf] = ComputePhaseM_SD_peaksSDF(MPeak0,SDPeak,TeoDur,condi,SO,ploing);
% 
% temdata = (Onset0.*360)./TeoDur; %Results in deegrees
% valsRad = deg2rad(temdata); %Results in radians
% [outMeanFF] = PhaseFanoFactor(valsRad,SALLFanoFactor,colorH,condi,SO,ploing);



ActivationRate.MRate = mean(SmoothRate);
ActivationRate.SDRate = std(SmoothRate);
ActivationRate.PeakRateMag = PeakRate;
ActivationRate.PeakRateOnset = Onset(BinPeakRate);%+TeoDur;
ActivationRate.PeakRateBin = BinPeakRate;
ActivationRate.Rate = SmoothRate;
ActivationRate.WeithedRateMean = wm;
ActivationRate.MiddleX = nnMb/2;
ActivationRate.Skewness = skeqq ;
ActivationRate.LG4MinAsyn = cf.A;
ActivationRate.LG4MaxAsyn = cf.D;
ActivationRate.LG4Slope = cf.B;
ActivationRate.LG4PSE = cf.C;
ActivationRate.LG4AdjR2 = G.adjrsquare;
ActivationRate.LG4Model = cf;
ActivationRate.NumThresh25 = numypred25;
ActivationRate.NumThresh75 = numypred75;
ActivationRate.NumThresh50 = numypred50;
ActivationRate.ypred = {ypred};

ActivationRate.totalCells = totalCells;


