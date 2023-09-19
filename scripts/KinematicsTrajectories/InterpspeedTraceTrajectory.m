function IntTrj = InterpspeedTraceTrajectory(Trj,NBinBase)
% NBinBase = 168;%168
[~, NBin] = size(Trj);
t = 1: NBin;
tBase = 1: (NBin-1)/(NBinBase-1): NBin;
IntTrj = spline(t, Trj(1, :), tBase);
IntTrj(IntTrj < 0 ) = 0; 