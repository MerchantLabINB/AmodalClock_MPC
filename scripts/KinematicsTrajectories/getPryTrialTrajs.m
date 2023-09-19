function Trajectories = getPryTrialTrajs(Manifold,MeanTrj)

V = Manifold.V_Dur;

for cEje = 1:3
    Resp = sum(MeanTrj(1:3, :)' .* repmat(V(:, cEje)', size(MeanTrj,2),1), 2);
    Trajectories(cEje,:) = Resp';
end