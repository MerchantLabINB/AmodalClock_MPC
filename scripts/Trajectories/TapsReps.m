function  TapsTrials = TapsReps(c_trials,SO)

%Estructura para Manifold
%---Matrix of 500 x 3-----%
%Taps x PC
Trajectories = c_trials;
Traj = cellfun(@(x) x(1:3,:),Trajectories,'un',0);%'
TapsTrials = [];
for i=1:length(c_trials)
    Binsize = size(Traj{i,1},2)/4;
    if(SO==5)
    movs = 0: Binsize : size(Traj{i,1},2);
    movs(1) = 1;    
    else
       movs = Binsize: Binsize : size(Traj{i,1},2);
%        movs(1) = movs(1)+1;
%        movs = movs(1:end-1);
    end
%     movs = movs  + Mov;  
%     TapsPapi = Traj{1,i}(:,movs - Mov2)';
    TapsPapi = Traj{i,1}(:,movs)';
    TapsTrials = [TapsTrials;TapsPapi];
    TapsPapi = [];
end