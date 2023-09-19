function oPeakQuarter = getPeakQuarterNum3(CellId,PeakQuarter,task,SO,TapTimes)

TapT = TapTimes(SO);
PeakQuarter = PeakQuarter-TapT;

if SO == 4
     TeoDur = [445 850 445 850];
else
    TeoDur = [445 850 445 850];
end

IntDur = TeoDur(task);

tQuarter = round(IntDur/4);
for ii = 1:4
   qs(ii,1) = (tQuarter*ii);
end
%qs(4,1) = qs(4,1)+60; %upper limit plus 60 of delta

Id1stQuarter =  CellId;
Id2ndQuarter =  CellId;
Id3rdQuarter =  CellId;
Id4rdQuarter =  CellId;
Id5thQuarter =  CellId;

Id1stQuarter(PeakQuarter > qs(1,1)) = [];
Id2ndQuarter(PeakQuarter <= qs(1,1) | PeakQuarter > qs(2,1)) = [];
Id3rdQuarter(PeakQuarter <= qs(2,1) | PeakQuarter > qs(3,1)) = [];
Id4rdQuarter(PeakQuarter <= qs(3,1) | PeakQuarter > 99800) = [];
Id5thQuarter(PeakQuarter <= 99800) = [];

oPeakQuarter(1,1) = numel(Id1stQuarter);
oPeakQuarter(1,2) = numel(Id2ndQuarter);
oPeakQuarter(1,3) = numel(Id3rdQuarter);
oPeakQuarter(1,4) = numel(Id4rdQuarter);
oPeakQuarter(1,5) = numel(Id5thQuarter);

