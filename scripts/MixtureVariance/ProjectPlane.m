function ProjectPlane(M,V,D, optPlane)

if(strcmp(optPlane,'Duration'))
    [u,v] = meshgrid(-1:2:1);
    %----Set Coordinates-----%+
    %Ancho
    v(1,1) = -1.45;  %1.7  %-1.4
    v(1,2) = -1.45; %1.7   %-1.4
    v(2,1) = 1.4;%1.2  %1.14
    v(2,2) = 1.4;%1.2
    
    % largo
    
    
    u(1,2) = 1.65;%1.8
    u(2,2) =1.65;  %1.8
    
    u(1,1) = -1.4;
    u(2,1) = -1.4;
    %Desviacion estandar
    for cDim = 1: 3
        X{cDim} = sqrt(D(2,2))*V(cDim, 2)*u + sqrt(D(3,3))*V(cDim, 3)*v;
    end
    % clf
    f1 = gcf
    set(0,'CurrentFigure',f1)
    hold on
    fPlane = 1.5;
    % Gráfica del plano ...
    % [0.94,0.97,1.00]
    mesh(fPlane*X{1}+M(1), fPlane*X{2}+M(2), fPlane*X{3}+M(3), zeros(size(X{3})), 'FaceColor',  [0.94,1,1.00],'FaceAlpha','0.15','LineWidth', 3)
else
    
    [u,v] = meshgrid(-1:2:1);%(-1:2:1)
    
    %-----Set Coordinates-----%
    
    u(1,2) = 1.02;%- 95  &1.18 Alto %1.5
    u(2,2) = 1.02;   %1.18      %Alto %1.5
    u(1,1) = -1.07;
    u(2,1) = -1.07;
    
    v(1,1) = -1;   %-1.4
    v(1,2) = -1;   %-1.4
    v(2,1) =1.03;%78  %1.14
    v(2,2) = 1.03;
    %desviacion estandar
    for cDim = 1: 3
        X{cDim} = sqrt(D(2,2))*V(cDim, 2)*u + sqrt(D(3,3))*V(cDim, 3)*v;
    end
    
    f1 = gcf;
    set(0,'CurrentFigure',f1)
    % clf
    hold on
    fPlane = 2;
    % Gráfica del plano ...
%     plot3(TrjCon(:, 1), TrjCon(:, 2), TrjCon(:, 3), '.k')
    mesh(fPlane*X{1}+M(1), fPlane*X{2}+M(2), fPlane*X{3}+M(3), zeros(size(X{3})), 'FaceColor',  [0.94,0.97,1.00],'FaceAlpha','0.15','LineWidth', 3)
    
end