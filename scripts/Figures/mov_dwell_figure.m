function mov_dwell_figure(mov_dwell_array,Pop)


close all
%-------Figure Parameters------%
f1 = figure;
whitebg('k')
set(f1, 'InvertHardcopy', 'off')
set(f1, 'color', [0 0 0]);

Cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);



%% Speed Profile
subplot(2,3,1:3)
plot(Pop(1, 4).Speed{1},'LineWidth',1.8)
xlabel('Time (ms)')
ylabel('Speed (a.u.)')
set(gca,'Ytick',[0:200:800])
set(gca,'Tickdir','out')
box off


%% Set movements and dwell
movements = mov_dwell_array.movs;
dwells = mov_dwell_array.dwells;
sizes = cellfun(@(x) size(x,1),movements,'un',1);
TaskLabel = {'A450','A850','V450','V850'};
mov_array = cat(1,movements{:});
dwell_array = cat(1,dwells{:});
grp=[]; %grouping matrix
dur_Array = []; %Duration array
short_large = [1,2,1,2];
for n=1:4
  grp=vertcat(grp,n*ones(sizes(n),1));%sizes is a variable with n column each column stores the size of each variable
  dur_Array = vertcat(dur_Array,short_large(n)*ones(sizes(n),1));
end
%% Create Subplot for
% Movement Profile and dwells


%% Box plot of movements
subplot(2,3,4)
ax = boxplot(mov_array,grp,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);
ApplyParams(gcf,Cmap)
ylabel('Movement duration (ms)')
ax = gca;
%% Box plot of dwells
subplot(2,3,6)
ax2 =boxplot(dwell_array,grp,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);
ApplyParams(gcf,Cmap)
ylabel('Dwell duration (ms) ')
ax2 = gca;

set([ax,ax2],'Ylim',[0 800])
set([ax,ax2],'YTick',0:200:800)

%% Anova
% 1 - Auditory 
% 2 - Visual
% 1 -  450
% 2 -  850
aud = ones(1,sum(sizes(1:2)));
vis =  ones(1,sum(sizes(3:4)))*2;
Modality = [aud,vis];

mov_stat = anovan(mov_array,{dur_Array';Modality},'model','full','varnames',{'Duration','Modality'},'display','on');
dwell_stat = anovan(dwell_array,{dur_Array';Modality},'model','full','varnames',{'Duration','Modality'},'display','on');
