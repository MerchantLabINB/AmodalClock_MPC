
function  Var_Axes_Trj= SerialOrder_Plane(Manifold,f)

G1 = [143, 101, 134]/255;
G2 = [200, 83, 82]/255;
G3 = [119, 38, 67]/255;
G4 = [71, 0, 30]/255;

cmapso = [G1;G2;G3;G4];

steps(1,:) = [1  42   84 126];
steps(2,:) = [42 84  126 168];



cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);

%We need to interpolate the trajectories for durations and conditions.
%We use same number of bins. We use 168 bins that corresond to the  number of bins for
% the large duration. Every 42 bins is when the mtapping occured.

Trj = cellfun(@InterpTrajectory, Manifold.trajectories, 'UniformOutput', 0);
% # Condition
% 1->A450, 2->A850, 3->V450, 4->V850
% for Cond = [4]
for Cond = 2
    
set(0,'CurrentFigure',f)

TrjCon = [Trj{:, Cond}]';
% Compute the Variability axes
[V, D] = eig(cov(TrjCon));
% La variabilidad explicada en estos ejes es
Var_Axes_Trj = diag(D) / sum(diag(D));
% El vector columna 1 de V es el vector tangente al plano, este vector es
% el que se usa para el c�lculo de �ngulo entre planos
% Los vectores columna 2 y 3 de V, determinan el plano
M = mean(TrjCon); % Promedio de datos
% C�lculo del plano
[u,v] = meshgrid(-1:2:1);
for cDim = 1: 3
    X{cDim} = sqrt(D(2,2))*V(cDim, 2)*u + sqrt(D(3,3))*V(cDim, 3)*v; 
end
% clf
% hold on
% plot3(TrjCon(:,1),TrjCon(:,2),TrjCon(:,3),'.','Color',[1 1 1],'MarkerSize',5)
hold on
for trial=1:25
for i=1:4
plot3(Trj{trial,Cond}(1,steps(1,i):steps(2,i)),Trj{trial,Cond}(2,steps(1,i):steps(2,i)),Trj{trial,Cond}(3,steps(1,i):steps(2,i)),'.','color',cmapso(i,:),'MarkerSize',10)
end
end



fPlane = 2.4;
% Gr�fica del plano ...
mesh(fPlane*X{1}+M(1), fPlane*X{2}+M(2), fPlane*X{3}+M(3), zeros(size(X{3})),'FaceColor',[.5,.5,.5],'EdgeColor',[0.94,0.97,1.00],'FaceAlpha','0.15', 'LineWidth', 3)
% los ejes ...
% ColorAxes = [1 0 0; 0 0 0; 0 0 0];
ColorAxes= repmat([0.94,0.97,1.00],3,1);
ColorAxes(1,:) = [1 1 0];
for cAxes = 1: 3
    line(M(1)+sqrt(D(cAxes,cAxes))*[0 V(1, cAxes)], M(2)+sqrt(D(cAxes,cAxes))*[0 V(2, cAxes)], M(3)+sqrt(D(cAxes,cAxes))*[0 V(3,cAxes)], 'color', ColorAxes(cAxes, :), 'linewidth', 2)
end
axis equal

% set(gca, 'CameraPosition', [7 11 50]);
% Variability(Cond,1) = Var_Axes_Trj(1) * 100;
% legend(['Variability - ',num2str(Variability(Cond,1))]);
% title(ConTitle{Cond})
% axis equal
% xlabel('PC 1','FontSize',12)
% ylabel('PC 2','FontSize',12)
% zlabel('PC 3','FontSize',12)
end
%-----GPFA-----%
% view([20,38])
%-----PCA V850-----%
% view([95 135 -100])
% view([135 -25])
% xlim([-150 100])
% zlim([-180 150])

