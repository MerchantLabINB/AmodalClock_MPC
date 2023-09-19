function ax = Conducta_Modelo_2021(Behave,f3)

%Set figure parameters
lwth = 3; %linewidth
ErrorStd = 5;%7
numplots = 4;
plots=3;%3

%From Behav Array
params.intervals{1} = [450,850,450,850];
MS1 = 10;
%----COLOR MAP-------%
CM = [0,0.7490,1; 0 0 1; 1 .5020 0;1 0 0];
t = 1;

C={Behave(:).ID}';
strTask = {'A450','A850','V450','V850'};
Ind = cell2mat(arrayfun(@(x) find(contains(C,strTask{x})), 1:4,'un',0));

bParams = cell(1);

SO = 1:5; %number of taps
endStim =5; %Find the second Stim

for k=1:4
    Taps = {Behave(Ind(:,k)).Movs};
    Taps = cellfun(@(x) x(1:25,SO),Taps,'un',0);
    
    Stims ={Behave(Ind(:,k)).Stims};
    Stims = cellfun(@(x) x(1:25,end-endStim:end-1),Stims,'un',0);%end-5 : end-1
    
    %-----Produced target interval-----%
    PItarget_Avg = cellfun(@(x) diff(x,[],2) *1000,Stims,'un',0);%ms 
    
    %-----Produced Interval-----%
    PITaps_Avg = cellfun(@(x) diff(x,[],2) *1000,Taps,'un',0);%ms
    
    PITrials = cellfun(@(x) mean(x,'omitnan'),PITaps_Avg,'un',0);
    PITrialsStd = cellfun(@(x) std(x,'omitnan'),PITaps_Avg,'un',0);
    PIAvgSO_STd = cat(1,PITrialsStd{:});
    PIAvgSO = cat(1,PITrials{:});
    
    %Lag -1
    PItarget_AvgSO =cellfun(@(x) mean(x),PItarget_Avg,'un',0);
    targetAvg = cat(1,PItarget_AvgSO{:});
    ConstantError = PIAvgSO - targetAvg;
    
 
    bParams{k,2} = PIAvgSO_STd(:);
    bParams{k,1} =  ConstantError(:);
    bParams{k,3} = targetAvg(:);
   
    %each row correspond with conditions
    %Auditory 450,850.  Visual 450,850.
    %22 Sessions
    %25 trials
    %4 Serial Order
    %Produced Interval = 25 * 4 * 22 = 2200
end


%-----------Plot Parameters----------%
%% Poner en funciones el codigo para graficar
set(0,'CurrentFigure',f3)

for j = 1:2
    ax(j) = subplot(plots,numplots,j);
    SEM = [];
    AvgTaskCond =[];
    
    for i = 1:4
        mPositions = bParams{i,j};
        
        if(j ==1)
            hold on;
            SEM(1,i) = std(bParams{i,1})/sqrt(length(bParams{i,1}));
            
            
            [ph1, ~]= boundedline(ax(j),params.intervals{t}(i),mean(mPositions),SEM(1,i),'o','cmap',CM(i,:),'alpha','transparency', 0.4);
            errorbar(params.intervals{t}(i),mean(mPositions),SEM(1,i)*ErrorStd,'color',CM(i,:),'linewidth',lwth);
        else
%             SEM(1,i) = std(bParams{i,2})/2;
              SEM(1,i) = std(bParams{i,2})/sqrt(length(bParams{i,2}));
            
            hold on;
            [ph1, ~]= boundedline(ax(j),params.intervals{t}(i),mean(mPositions),SEM(1,i),'o','cmap',CM(i,:),'alpha','transparency', 0.4);
            errorbar(params.intervals{t}(i),mean(mPositions),SEM(1,i)*ErrorStd,'color',CM(i,:),'linewidth',lwth);
            
        end
        set(ph1,'MarkerSize',MS1,'MarkerEdgeColor',CM(i,:),'MarkerFaceColor',CM(i,:));
    end   
end





clear SEM
[LagOutput, LagTrials] = autocorr_lag1(Behave);
Lags.LagOutput = LagOutput;
Lags.LagTrials = LagTrials;

j = 4;
ax(j) = subplot(plots,numplots,j);

for i = 1:4
    hold on
    mPositions = LagOutput(:,i);
    SEM(1,i) = std(mPositions)/sqrt(length(mPositions));
    [ph1, ~]= boundedline(ax(j),params.intervals{t}(i),mean(mPositions),SEM(1,i),'o','cmap',CM(i,:),'alpha','transparency', 0.4);
    set(ph1,'MarkerSize',MS1,'MarkerEdgeColor',CM(i,:),'MarkerFaceColor',CM(i,:));
    errorbar(params.intervals{t}(i),mean(mPositions),SEM(1,i)*ErrorStd,'color',CM(i,:),'linewidth',lwth);
    
end
ax(3) = subplot(plots,numplots,3);
plot(nan,nan)

hold on
x=1:900;
y=0;
plot(ax(1),x,y*ones(size(x)),'Color','w','LineStyle','--')
plot(ax(4),x,y*ones(size(x)),'Color','w','LineStyle','--')

