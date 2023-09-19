function  [Vp,VariabilityAxes,VarAcum] = subspacePlane(TrjCon,opt,optPlane)


[V, D] = eig(cov(TrjCon));
% V = V(:, [3, 2, 1]);
Vp(1,:) = V(:,1); %Eje del plano tangente
% La variabilidad explicada en estos ejes es
Var_Axes_Trj = diag(D) / sum(diag(D));
VariabilityAxes = Var_Axes_Trj;
VarAcum = sum(VariabilityAxes(2:3));
M = mean(TrjCon); % Promedio de datos
% ColorAxes = [1 0 0; .5 .5 .5; .5 .5 .5];





if(opt == 1)
%     clear Vp
    Vp = V;
    ProjectPlane(M,V,D, optPlane)

% ColorAxes = repmat([255 240 31]/255,3,1);
% ColorAxes= repmat([0.94,0.97,1.00],3,1);
% figure;
% for cAxes =  3
%     line(M(1)+sqrt(D(cAxes,cAxes))*[0 V(1, cAxes)], M(2)+sqrt(D(cAxes,cAxes))*[0 V(2, cAxes)], M(3)+sqrt(D(cAxes,cAxes))*[0 V(3,cAxes)], 'color', ColorAxes(cAxes, :), 'linewidth', 2)
% 
% end
end