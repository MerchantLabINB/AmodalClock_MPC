function [portionArea_dwell,movs_segment]= getAUC_movs(params,binsdwell,cModality,sp)
% addShift = [28;28;32;29]-1;
step = sp(cModality);
sp2 = [0,step:step:step*3];

%Restamos el periodo control
Distances = params{cModality};
xx = binsdwell; %26

if(xx(1,1) <= 0 )
    xx(1,1)= 1;
end
for numTrials = 1:size(params{cModality},1)
    aux_Trj = [];
    traj = Distances(numTrials,:);

for j = 1:4

    reTraj = reshape(traj(~isnan(traj)),sp(cModality),4);
x= 1:length(traj); 
y = traj;


A = cumtrapz(x, y);
newFunction = @(p,q) max(A(x<=q)) - min(A(x>=p));
portionArea_dwell(numTrials,j) = newFunction(xx(j,1), xx(j,2));

Trj = y(xx(j,1):xx(j,2));
IntTrj = InterpspeedTraceTrajectory(Trj,42);
aux_Trj = [aux_Trj,IntTrj];

end
movs_segment(numTrials,:)  =aux_Trj;


end

