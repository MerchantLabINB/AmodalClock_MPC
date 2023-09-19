function PeakSkew = ComputePeakSkewAPDur(Peak,Begin,Dur)


DifPeakBeg = Peak -Begin;
PeakSkew = DifPeakBeg./Dur;