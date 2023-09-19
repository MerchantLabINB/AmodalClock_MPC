function [SortedCells,CellGroups,PeakQuarterNum] = getSortedCellsALLSO(Tpeak,TapTimes,task,Serial_order)

ThresholdSD = 360;%160
MinimumSOAP = 2;
ThresholdMedian = 150;
TapTimes = round(TapTimes);

CellGroups = ones(numel(Tpeak(:,1)),1);

TotAPSO(:,1) = 4 - sum(isnan(Tpeak'));

PeakAligned = Tpeak - TapTimes(1:4);
OriginalCells(:,1) = 1:numel(Tpeak(:,1));

tPeakAligned = PeakAligned;
tPeakAligned(TotAPSO<MinimumSOAP,:) = 99999;
CellGroups(TotAPSO<MinimumSOAP,:) = 3; %%only one activation period

MedianPeak(:,1) = nanmedian(tPeakAligned');
ttPeakAligned = tPeakAligned;

NumNans = zeros(numel(Tpeak(:,1)),1);
for so =1:4    
  difmedian =  abs(ttPeakAligned(:,so)- MedianPeak);
  ttPeakAligned(difmedian > ThresholdMedian,so) = NaN;
  ttPeakAligned(TotAPSO<3,so) = tPeakAligned(TotAPSO<3,so);
  NumNans = NumNans+ (difmedian > ThresholdMedian);
end
ttPeakAligned(NumNans == 4,:) =  99900;%%%cells with high variability
CellGroups(NumNans == 4,:) =  2;%%%cells with high variability

MedianPeak(:,1) = nanmedian(ttPeakAligned');
SDPeak(:,1) = nanstd(ttPeakAligned');


MedianPeak(SDPeak > ThresholdSD) = 99900;
CellGroups(SDPeak > ThresholdSD) = 2;   %%%high variability

[SortMedian,SortedCells] =  sort(MedianPeak);

sorT4 = ttPeakAligned(SortedCells,:);

temcelgro = CellGroups(SortedCells);

  tCellIds = 1:1:numel(SortMedian);
  PeakQuarterNum = getPeakQuarterNum2(tCellIds,SortMedian,task,Serial_order);







