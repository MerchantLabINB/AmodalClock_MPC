function ax = mov_dwell_trace_panels(SpeedParams,f1,ax)

mov_dwell_array.movs = cellfun(@(x) x(:),SpeedParams.mov_SO,'un',0);
mov_dwell_array.dwells = cellfun(@(x) x(:),SpeedParams.dwell_SO,'un',0);


SpeedParams.mean_trace_mov

Cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);

set(0,'CurrentFigure',f1)

cols = 5;
rows = 8;

%% Speed Profile Panel H
% figure;
ax(10) = subplot(rows,cols,[14,19]);

hold on;
for i=1:4
Speed_new = SpeedParams.mean_trace_mov{1,i};
plot(Speed_new,'LineWidth',1.8,'Color',Cmap(i,:));
end

%Panel I
ax(11) = subplot(rows,cols,[15,20]);
hold on;
for i=1:4
plot(SpeedParams.mean_trace_dwell{1,i},'LineWidth',1.8,'Color',Cmap(i,:))
end

Tapp = [25 , 25];
ymax = [0 70];%get(ax_mv(1),'YLim');
plot(ax(10),[Tapp; Tapp],[ymax(1) ymax(2)],'--','Color',[0 , 1, 0],'LineWidth',1.4)

%% Set movements and dwell
movements = mov_dwell_array.movs;
dwells = mov_dwell_array.dwells;
sizes = cellfun(@(x) size(x,1),movements,'un',1);
sizes2 = cellfun(@(x) size(x,1),dwells,'un',1);

mov_array = cat(1,movements{:});
dwell_array = cat(1,dwells{:});
grp=[]; %grouping matrix
grp2=[]; %grouping matrix

dur_Array = []; %Duration array
dwell_Array = []; %Dwell array
short_large = [1,2,1,2];
for n=1:4
  grp=vertcat(grp,n*ones(sizes(n),1));%sizes is a variable with n column each column stores the size of each variable
  grp2=vertcat(grp2,n*ones(sizes2(n),1));%sizes is a variable with n column each column stores the size of each variable
  dur_Array = vertcat(dur_Array,short_large(n)*ones(sizes(n),1));
  dwell_Array = vertcat(dwell_Array,short_large(n)*ones(sizes2(n),1));
end

%% Box plot of movements
ax(12) = subplot(rows,cols,[26,31]);
TaskLabel = {'AS','AL','VS','AL'};
boxplot(ax(12),mov_array,grp,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);
ApplyParams(ax(12),Cmap)
%% Box plot of dwells
ax(13) = subplot(rows,cols,[27,32]);
TaskLabel = {'AS','AL','VS','VL'};
boxplot(ax(13),dwell_array,grp2,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);
ApplyParams(ax(13),Cmap)
%Population

%% Anova
% 1 - Auditory 
% 2 - Visual
% 1 -  450
% 2 -  850
% aud = ones(1,sum(sizes(1:2)));
% vis =  ones(1,sum(sizes(3:4)))*2;
% Modality = [aud,vis];
% 
% aud = ones(1,sum(sizes2(1:2)));
% vis =  ones(1,sum(sizes2(3:4)))*2;
% Modality2 = [aud,vis];
% 
% % mov_stat = anovan(mov_array,{dur_Array';Modality},'model','full','varnames',{'Duration','Modality'},'display','off');
% % dwell_stat = anovan(dwell_array,{dwell_Array';Modality2},'model','full','varnames',{'Duration','Modality'},'display','on');
% this_stat_dwell = KinematicsStatistics(dwells);
% this_stat_mov = KinematicsStatistics(movements);
% 
% % %% Articulo ANOVA
% Speed_dwell = cellfun(@(x) x(:,2:end),SpeedParams.dwell_SO,'UniformOutput',0);
% Speed_mov = cellfun(@(x) x(:,2:end-1),SpeedParams.mov_SO,'UniformOutput',0);
% this_stat_dwell = KinematicsStatistics(Speed_dwell);
% this_stat_mov = KinematicsStatistics(Speed_mov);
% 
% 
% 
% this_stat = KinematicsStatistics(Speed_mov,Speed_dwell);