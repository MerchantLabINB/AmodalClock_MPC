function ax = mov_dwell_figure(mov_dwell_array,Pop,f3,ax)


Cmap = ([0,0.7490,1;...
    0 0 1;...
    1 .5020 0;...
    1 0 0]);

    set(0,'CurrentFigure',f3)

%% Speed Profile
ax(5) = subplot(3,4,6:7);
Taps = Pop(1, 1).TimesTap(1,2:7);%32
VidSpeed = Pop(1, 1).Speed{1};%32
TimesVideo = Pop(1, 1).Times{1};
% VidSpeed = VidSpeed(TimesVideo >= Taps(1)-.06 & TimesVideo <= Taps(end)+.15);
% VidTimes = TimesVideo(TimesVideo >= Taps(1)-.06 & TimesVideo <= Taps(end)+.15);
VidSpeed = VidSpeed(TimesVideo >= Taps(1)-.2 & TimesVideo <= Taps(end)+.2);
VidTimes = TimesVideo(TimesVideo >= Taps(1)-.2 & TimesVideo <= Taps(end)+.2);
plot(VidTimes*1000,VidSpeed)
% box off
%% Set movements and dwell
movements = mov_dwell_array.movs;
dwells = mov_dwell_array.dwells;
sizes = cellfun(@(x) size(x,1),movements,'un',1);
sizes2 = cellfun(@(x) size(x,1),dwells,'un',1);

TaskLabel = {'AS','AS','VS','VL'};
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
%% Create Subplot for
% Movement Profile and dwells


%% Box plot of movements
ax(6) = subplot(3,4,10);
boxplot(ax(6),mov_array,grp,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);

ApplyParams(ax(6),Cmap)

% ax = gca;
%% Box plot of dwells
ax(7) = subplot(3,4,11);
boxplot(ax(7), dwell_array,grp2,'notch','on','colors',Cmap,'Whisker',0,'Symbol','',...
'labels',TaskLabel);
ApplyParams(ax(7),Cmap)
% ApplyParams(gcf,Cmap)




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
% mov_stat = anovan(mov_array,{dur_Array';Modality},'model','full','varnames',{'Duration','Modality'},'display','off');
% dwell_stat = anovan(dwell_array,{dwell_Array';Modality2},'model','full','varnames',{'Duration','Modality'},'display','off');
