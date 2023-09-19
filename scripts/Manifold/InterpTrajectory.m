
function IntTrj = InterpTrajectory(Trj)
NBinBase = 168;%168;%235;
[~, NBin] = size(Trj);
t = 1: NBin;
tBase = 1: (NBin-1)/(NBinBase-1): NBin;
IntTrj = spline(t, Trj(1:3, :), tBase);
% IntTrj = spline(t, Trj(1:1019, :), tBase);
% IntTrj = spline(t, Trj(1:end, :), tBase);
