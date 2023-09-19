    function [Trajectories,TapsPry,TapsRespPry] = getProyections(Manifold)
    
    
    V = Manifold.V_Dur;
    
    % TrajectoriesDur = Manifold.trajectories;
    AllTaps = Manifold.TapsReps;
    len = size(AllTaps,1);
    TapsPry(:,1) = sum(AllTaps.* repmat(V(:, 1)', len,1),2);
    TapsPry(:,2) = sum(AllTaps.* repmat(V(:, 2)', len,1),2);
    TapsPry(:,3) = sum(AllTaps.* repmat(V(:, 3)', len,1),2);
    
    MeanTrj = cellfun(@InterpTrajectory, Manifold.Meantrajectories, 'UniformOutput', 0);
    Movs = 0:42:168;
    
    Movs(1) = 1;
    for cEje = 1:3
        for cMod = 1: 4
                Resp = sum(MeanTrj{cMod}(1:3, :)' .* repmat(V(:, cEje)', 168,1), 2);
               TapsPryPoint = MeanTrj{cMod}(1:3, :)' .* repmat(V(:, cEje)', 168,1);
    
            TapsRespPry{cMod}= TapsPryPoint(Movs',:);
            Trajectories{cMod}(cEje,:) = Resp';
        end
        
    end
    