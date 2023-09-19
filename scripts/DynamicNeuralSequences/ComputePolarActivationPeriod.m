function OutCircStatistics = ComputePolarActivationPeriod(Begin,Dur,TeoDur,condi,SO,TapTime,ploing,PlotOrder)
%%Onset of Peak

Begin(Begin == 99999) = [];
Dur(Dur== 99999) = [];

Tneurons = numel(Begin);
Begin = round(Begin);
Dur = round(Dur);
TapTime = round(TapTime);
colorH = 2;

Begin = Begin - TapTime;
i = 1;
for neu = 1:Tneurons
    bg = Begin(neu);
    durt = Dur(neu);
    for ii = bg:bg+durt
       ActivationTime(i,1) = ii;
       i = i+1;
    end
    MeanAP(neu,1) = round(mean(ActivationTime));
end


%[OutCircStatistics] = ComputePeakPolarDistri(ActivationTime,TeoDur,2*pi,2*pi,condi,SO,ploing,PlotOrder,0);
[OutCircStatistics] = ComputePeakPolarDistri(ActivationTime,TeoDur,2*pi,2*pi,condi,ploing,PlotOrder);
                                                                      
% MeanR_degrees = rad2deg(r_angle);
% if MeanR_degrees < 0
%     MeanR_degrees = 360+MeanR_degrees ;
% end




% PolarActivation.MeanResultant = MeanResultant;
% PolarActivation.MeanR_degrees = MeanR_degrees;
% PolarActivation.pvalRalegigh = pvalRalegigh;
% PolarActivation.outMeanR  = outMeanR;



