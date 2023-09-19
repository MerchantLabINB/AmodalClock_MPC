function [portionArea_dwell,TotalArea_dwell,PeaktoPeak_Amp,dwell_segment]= getAUC(params,binsdwell,cModality,sp)

step = sp(cModality);
sp2 = [0,step:step:step*3];

%Restamos el periodo control
Distances = params{cModality};
xx = binsdwell;

if(xx(1,1) < 0 )
    xx(1,1)= 1;
end

for numTrials = 1:size(params{cModality},1)
    aux_Trj = [];
    traj = Distances(numTrials,:);
    x= 1:length(traj);
    y = traj;

    reTraj = reshape(traj(~isnan(traj)),sp(cModality),4);

    for j = 1:4
        ySo = reTraj(:,j)';
        A = cumtrapz(x,y);
        newFunction = @(p,q) max(A(x<=q)) - min(A(x>=p));
        try


            portionArea_dwell(numTrials,j) =newFunction(xx(j,1),xx(j,2));
            TotalArea_dwell(numTrials,j) = trapz(1:size(ySo,2),ySo);
            PeaktoPeak_Amp(numTrials,j) = peak2peak(y);

            Trj = y(xx(j,1):xx(j,2));
            IntTrj = InterpspeedTraceTrajectory(Trj,42);
            aux_Trj = [aux_Trj,IntTrj];

        end
    end
    dwell_segment(numTrials,:)  =aux_Trj;
end

