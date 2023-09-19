function [Trajectories,TapsPry,TapsRespPry] = getProyectionsVariable(Manifold)


V = Manifold.V_Dur;
AllTaps = Manifold.TapsReps;
len = size(AllTaps,1);

TapsPry(:,1) = sum(AllTaps.* repmat(V(:, 1)', len,1),2);
TapsPry(:,2) = sum(AllTaps.* repmat(V(:, 2)', len,1),2);
TapsPry(:,3) = sum(AllTaps.* repmat(V(:, 3)', len,1),2);


MeanTrj = Manifold.Meantrajectories;
Movs(1,:) = [1,22:22:88];
Movs(2,:) = [1,42:42:168];
Movs = repmat(Movs,2,1);

for cEje = 1:3

    for cMod = 1: 4
        Resp = sum(MeanTrj{cMod}(1:3, :)' .* repmat(V(:, cEje)', size(MeanTrj{cMod},2),1), 2);
        TapsPryPoint = MeanTrj{cMod}(1:3, :)' .* repmat(V(:, cEje)',size(MeanTrj{cMod},2),1);
        if(cEje == 1)
        end
        TapsRespPry{cMod}= TapsPryPoint(Movs(cMod,:)',:);
        Trajectories{cMod}(cEje,:) = Resp';
    end

end














