function [Vp,VariabilityAxes,XYZ] = ModalityPlane(f1,Manifold)
% Para poder comparar la información entre duraciones,
% vamos a interpolar las trayectorias con el mismo número de subintervalos,
% tomaremos 168 como base, que corresponde al número de subintervalos para
% la duración larga. Por lo tanto los momentos del tap, ocurriran ambos en
% cada 42 subintervalos
Trj = cellfun(@InterpTrajectory, Manifold.trajectories, 'UniformOutput', 0);
% Número de Condición
% 1->A450, 2->A850, 3->V450, 4->V850
Planes = [1,2;3,4];
% Planes = [1,2,3,4];
XYZ = [];
for aux = 1:2
% for aux = 1

Cond = Planes(aux,:);
TrjCon = [Trj{:, Cond}]';
% Cálculo de los ejes de variabilidad
[V, D] = eig(cov(TrjCon));
Vp(aux,:) = V(:,1);
% La variabilidad explicada en estos ejes es
Var_Axes_Trj = diag(D) / sum(diag(D));
% El vector columna 1 de V es el vector tangente al plano, este vector es
% el que se usa para el cálculo de ángulo entre planos
% Los vectores columna 2 y 3 de V, determinan el plano
M = mean(TrjCon); % Promedio de datos
% Cálculo del plano
% figure;
[u,v] = meshgrid(-1:2:1);%(-1:2:1)

%-----Set Coordinates-----%

u(1,2) = 1.02;%- 95  &1.18 Alto %1.5
u(2,2) =1.02;   %1.18      %Alto %1.5
u(1,1) = -1.07;
u(2,1) = -1.07;

v(1,1) = -1;   %-1.4
v(1,2) = -1;   %-1.4
v(2,1) =1.03;%78  %1.14
v(2,2) = 1.03;  
%------GPFA-----%
% -------Alto------%
% u(1,2) = 8.;%- 95  &1.18 Alto %1.5
% u(2,2) =.85;   %1.18      %Alto %1.5
% u(1,1) = -1.2;
% u(2,1) = -1.2;
% %-----Ancho------%
% v(1,1) = -1.1;   %-1.4
% v(1,2) = -1.1;   %-1.4
% v(2,1) = .95;%78  %1.14
% v(2,2) = .95;    
% 
% 


% if(aux == 1)
% u = u*1.30;
% u(1,1) = -1.4;%- 95
% % u(2,1) =-1.4;
% v = v* 1.08; %1.04
% else
% u = u*1.2;%1.09 %Ancho del rectangulo
% v(1,1) = -1.18;
% v(1,2) = -1.18;
% v(2,1) = .92;%78
% v(2,2) = .92;
% end
%desviacion estandar
for cDim = 1: 3
    X{cDim} = sqrt(D(2,2))*V(cDim, 2)*u + sqrt(D(3,3))*V(cDim, 3)*v; 
end

set(0,'CurrentFigure',f1)
% clf
hold on
fPlane = 2;
% Gráfica del plano ...
mesh(fPlane*X{1}+M(1), fPlane*X{2}+M(2), fPlane*X{3}+M(3), zeros(size(X{3})), 'FaceColor',  [0.94,0.97,1.00],'FaceAlpha','0.15','LineWidth', 3)
% los ejes ...
% ColorAxes = [1 0 0; .5 .5 .5; .5 .5 .5];
% ColorAxes= repmat([0.94,0.97,1.00],3,1);
% 
% for cAxes = 1: 3
%     line(M(1)+sqrt(D(cAxes,cAxes))*[0 V(1, cAxes)], M(2)+sqrt(D(cAxes,cAxes))*[0 V(2, cAxes)], M(3)+sqrt(D(cAxes,cAxes))*[0 V(3,cAxes)], 'color', ColorAxes(cAxes, :), 'linewidth', 2)
% end
xx = fPlane*X{1}+M(1);
yy = fPlane*X{2}+M(2);
zz = fPlane*X{3}+M(3);
XYZ = [XYZ; xx(:),yy(:),zz(:)];
VariabilityAxes(:,aux) = Var_Axes_Trj;
end
% view([20 -12])%PCA

% Modalidad  20 -12