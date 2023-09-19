function  plotRoutes3D_Updated(params,routes, movs,h1)
% Binslag = repmat([6 3],1,2);
Binslag = repmat([0 0],1,2);

binwidths = routes;
whitebg('k'); %Comentar para ver si se puede coparr figura

cmap=colormap();
% close
independent_windows=1;

Plot.Taps = 1;
Plot.headers = 1;
Plot.Points = 0;
Plot.Phase = params.Phase; 
Manifold = params.Manifold;


if(independent_windows)
%     h1=figure;
%     figure(h1);
    set(h1, 'InvertHardcopy', 'off')
    set(gcf, 'color', [0 0 0]);
    
else
    subplot(12,1,1:9)
end

% set(gca,'Linewidth',2);
% set(gca,'TickDir','out');

if(exist('binwidths'))
    %arrTimes=cellfun(@mean,binwidths,'UniformOutput',false); %aqui
    
    arrTimes=binwidths;%aqui para no simult
else
    arrTimes=[];
end

%%%%%Here you specify the width of he trajectory
%and markers.

%----PCA Firing Rate-----%
% StartEnd = .2;%1.5 3.5
% headerSize =.2;%2Sinlge session

% linewd =1.8;%PCA-.05 1.35
% linewd = 1;%PCA-.05 1.35
% linewd = 1.5;%CAPI
%-----LabSort------%
% linewd = 10.3;% 
% linewd = 0.02;%%Oswie coeff = 0.2
%------Kilosort-------%
% StartEnd = 6;%1.5 3.5
% linewd = 100;%CAPI
% headerSize = 6;


% %-----PCA-LabSort SubPopulation----%
% StartEnd = 1;%1.5 3.5
% linewd = 1;%CAPI
% headerSize = 1;
% facAmp = 4;

%Subpop withou cells
StartEnd = .5;%1.5 3.5
linewd = .5;%CAPI
headerSize = 1;
facAmp = 3;


%-----PCA-LabSort Population----%
% StartEnd = .1;%1.5 3.5
% linewd = .05;%CAPI
% headerSize = .1;
% facAmp = 1;
% facAmp = 2;
%-----GPFA----%
% % StartEnd = 3;%1.5 3.5
% % linewd = 3;%CAPI
% % headerSize = 3;

plotFunc=true;

ftimeSteps=cellfun(@arrTimeSteps,arrTimes,'UniformOutput',false);

if(exist('limTime'))
    maxBins=cellfun(@findTimeBin,ftimeSteps, num2cell(repmat(limTime,1,size(ftimeSteps,2))),'UniformOutput',false);
    if(isempty(maxBins{1}))
        plotFunc=0;
    else
        if(min(cell2mat(maxBins))<0)
            plotFunc=0;
        else
            routes=cellfun(@cutTrajectory,routes,num2cell(maxBins), 'UniformOutput',false);
        end
    end
else
    maxBins=num2cell(cellfun(@(x) size(x,2),ftimeSteps));
    
end

% PLOT ONLY Synch OR CONT
minBins={1 1 1 1 };

%NORMAL
range=(1.-1./logspace(0,1,250)+.1);
blueScheme=horzcat(zeros(250,1),zeros(250,1),range');
redScheme=horzcat(range',zeros(250,1),zeros(250,1));
dodgerBlueScheme = repmat([0,0.7490,1],250,1);
orangeredScheme = repmat([1 .5020 0],250,1);


% cmapoints = [0,0.7490,1; 0 0 1; 1 .5020 0;1 0 0];

%Range of Colors by Condition
%15 250;     %A450
%546 680;    %A850
%765 1000;   %V450
%1015 1250   %V850
%----------------------------%
task = params.tasks{1};
%Volver a descomentar
if(task(1)=='m')
    colorLims=[15 250;260 500;515 750;765 1000;1015 1250;1265,1500]; %MTV
    
else
    colorLims=[15 250;546 680;765 1000;1015 1250]; %Blocks Task
end

%   opengl hardware <-Mejor visualizacion
set(gca,'Color',[0 0 0]); %#2

% idx = [1];
if(plotFunc)
    for idx = 1:4%length(intervals) %aqui 1:4:length para saltar intervalos
        if(independent_windows)
%             hMain=figure(h1);
            hMain = h1;
        else
            hMain=subplot(12,1,1:9);
        end
        set(gca,'Color',[0 0 0]);%Descomentar
%         set(gca,'Linewidth',2);
%         set(gca,'TickDir','out');
        try
            p=routes{idx}(:,minBins{idx}:maxBins{idx});
        catch me
            continue;
        end
        
        
        curmov=ceil(movs{idx}); %aqui para no simutan
        if(exist('limTime'))
            curstim=cutStims(curstim,maxBins(idx));
        end
        
        totTS=size(p,2);
        v=axis();
        
        colorLine=linspace(colorLims(idx,1),colorLims(idx,2),totTS);
        
        %mainp=clinep(p(1,:),p(2,:),p(3,:),colorLine(1:totTS),linewd);
        %                 whitebg('k');
        set(gca, 'XGrid', 'off');
        set(gca, 'YGrid', 'off');
        set(gca, 'ZGrid', 'off');
        daspect([1 1 1]);
        [sx,sy,sz]=tubeplot(p(params.dims(1),:),p(params.dims(2),:),p(params.dims(3),:),linewd,colorLine(1:totTS),100);%,10-PCA,[1 0 1]); %[1 0 1]
        mainp=surf(sx,sy,sz,repmat(colorLine(1:totTS)',1,size(sx,2)));
        set(mainp,'EdgeColor','none','FaceLighting','phong', ...
            'EdgeLighting','phong','DiffuseStrength',1, ...
            'SpecularStrength',.11,'AmbientStrength',1);
        hold on;
        if (Plot.Points)
            %              metallicgold = [ 212 175 55 ] /255;
            %             amarant = [229,43,80]/255;
            %         	 amarant = [192, 43, 229] / 255;
            %              goldenyellow = [ 255 223 0 ] /255;
            %         	 amarant = [253, 204, 13] / 255;
            amarant = [255, 255, 0] / 255;
            
            %         %------Point in State Phase-----%
            dp = mean(cell2mat(SOMean'));
            [V,F]=platonic_solid(1,linewd*9);%Cube 2.8 Abraham
            orV=repmat(dp,size(V,1),1);%1
            patch('Faces',F, ...
                'Vertices',V+orV, ...
                'FaceVertexCData',repmat(amarant,size(F,1),1),'FaceColor','flat','FaceAlpha',0.9,'FaceLighting','phong','EdgeLighting','phong');
            
            %-----Starting Point Control-----%
            % %             SpC= SPControl{1,idx}(1:3,1)';
            % %             [V,F]=platonic_solid(2,linewd*4.5);%Cube 2.8 Abraham
            % %             orV=repmat(SpC,size(V,1),1);%1
            % %             patch('Faces',F, ...
            % %                 'Vertices',V+orV, ...
            % %                 'FaceVertexCData',repmat(cmapoints(idx,:),size(F,1),1),'FaceColor','flat','FaceAlpha',0.9,'FaceLighting','phong','EdgeLighting','phong');
        end
        %
        %
       if(Plot.Phase)
            %Head
%             greenTap = [0 255 0]/255; %PCA Size= 3.5
            goldenyellow = [0.5 0.5 0.5] ;

            [V,F]=platonic_solid(5,headerSize*StartEnd+1.2);%Cube 2.8 Abraham
            orV=repmat([p(params.dims(1),Plot.Phase(idx)),p(params.dims(2),Plot.Phase(idx)),p(params.dims(3),Plot.Phase(idx))],size(V,1),1);%1
            patch('Faces',F, ...
                'Vertices',V+orV, ...
                'FaceVertexCData',repmat(goldenyellow,size(F,1),1),'FaceColor','flat','FaceAlpha',0.9,'FaceLighting','phong','EdgeLighting','phong');           
       end
        
        if(Plot.headers)
            %Head
            greenTap = [250 250 210]/255; %PCA Size= 3.5
            [V,F]=platonic_solid(2,headerSize*StartEnd);%Cube 2.8 Abraham
            orV=repmat([p(params.dims(1),1),p(params.dims(2),1),p(params.dims(3),1)],size(V,1),1);%1
            patch('Faces',F, ...
                'Vertices',V+orV, ...
                'FaceVertexCData',repmat(greenTap,size(F,1),1),'FaceColor','flat','FaceAlpha',0.9,'FaceLighting','phong','EdgeLighting','phong');
            
            %End
            [V,F]=platonic_solid(3,headerSize*StartEnd); %Octahedron
            orV=repmat([p(params.dims(1),end),p(params.dims(2),end),p(params.dims(3),end)],size(V,1),1); %totTs
            patch('Faces',F, ...
                'Vertices',V+orV, ...
                'FaceVertexCData',repmat(greenTap,size(F,1),1),'FaceColor','flat','FaceAlpha',0.9,'FaceLighting','phong','EdgeLighting','phong');
        end
        
        k=6;
        n=2^k-1;
        hold on;
        
        bulColor=repmat([0 255 0]/255,6,1);
        plotMov = curmov; %Continous Time
        plotMov(2:end) = curmov(2:end)-Binslag(idx);
        plotMov = plotMov(2:end);
        %Aqui grafico mis taps
        if(Plot.Taps)
            for plotCMov=1:length(plotMov)
                
                [x, y, z] = ellipsoid(p(params.dims(1),plotMov(plotCMov)), ...
                    p(params.dims(2),plotMov(plotCMov)), ...
                    p(params.dims(3),plotMov(plotCMov)), ...
                    linewd*facAmp,linewd*facAmp,linewd*facAmp,n);
%                     linewd*2,linewd*2,linewd*3.2,n);

                for x1=1:size(x,1)
                    for x2=1:size(x,1)
                        cColor(x1,x2,1:3)=bulColor(plotCMov,:);
                    end
                end
                ps=surf(x, y, z,cColor);
                
                set(ps,'FaceAlpha',0.9,'EdgeColor','none');
                
            end
        end
        %Change 3D point of view
        if(exist('viewcoords'))
            campos(viewcoords);
            %Ajustar para cada pelicula
            axis(params.plotAxis);
        end
%         xlabel(['PC ' num2str(params.dims(1))],'FontSize',12);
%         ylabel(['PC ' num2str(params.dims(2))],'FontSize',12);
%         zlabel(['PC ' num2str(params.dims(3))],'FontSize',12);
%         
        cmap=zeros(1250,3);
        cmap(1:250,:) = dodgerBlueScheme(1:250,:);    %A450
        cmap(501:750,:)=blueScheme(1:250,:);          %A850
        cmap(751:1000,:) = orangeredScheme(1:250,:);  %V450
        cmap(1001:1250,:)=redScheme(1:250,:);         %V850
        cmap(1,:)=[1 1 1];
        
        colormap(cmap);
        
    end
    
    if(independent_windows)
%         figure(h1);
    else
        subplot(12,1,1:9);
    end
    a=axis();
    numSteps=length(routes{1});
    clinep(1:1250,zeros(12500,1)+1000000,zeros(1250,1)+1000000,1:1250);
    axis(a)
    
    material 'metal';
    light('Position',[0 0 1],'Style','infinite');
    light('Position',[0 0 -1],'Style','infinite');
    light('Position',[0 1 0],'Style','infinite');
    light('Position',[0 -1 0],'Style','infinite'); %Esta vista
    light('Position',[6 0 -1],'Style','infinite');
    
    
%     colormap(cmap);
%     if(independent_windows)
%         c=caxis();
%         caxis([0,c(end)]);
%     else
%         caxis(caxis(hMain));
%     end
    %-----MANIFOLD-----% 
    %absolet
%     if(params.Plane == [0 0 1] )%DURATION
%         [Vp,VariabilityAxes] = Duration_Plane(h1,Manifold);
%         XYZ = [];
%     elseif(params.Plane == [0 1 0] )%MODALITY
%         [Vp ,VariabilityAxes,XYZ] = ModalityPlane(h1,Manifold);
%     else
%       Vp  = [];
%       VariabilityAxes = [];
%       XYZ = [];
%     end
    %-------------------%    
end
% OutParams.XYZ = XYZ;
% OutParams.Variability = VariabilityAxes;
% OutParams.Vp = Vp;
end


function timesteps=arrTimeSteps(arrTimes)
%    lastTime=0;
%    for cont=1:length(arrTimes)
%        timesteps(cont)=arrTimes(cont)+lastTime;
%        lastTime=timesteps(cont);
%    end
timesteps=cumsum(arrTimes);
end

function bin=findTimeBin(timesteps, time)
bin= max(find(timesteps<time))
if(isempty(bin) | bin<2 )
    bin=2;
end
end

function newTrajectory=cutTrajectory(orgTrajectory,maxBin)
if(iscell(maxBin))
    maxBin=maxBin{1};
end
newTrajectory=orgTrajectory(:,1:maxBin);
end

function newStims=cutStims(orgStims,maxBin)
if(iscell(maxBin))
    maxBin=maxBin{1};
end
newStims=orgStims(1:max(find(orgStims<maxBin)));
end