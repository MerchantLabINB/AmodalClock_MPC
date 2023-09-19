function outputStruct = get_trajectoriesPopulation(outputStruct)
% clear all
params.projection=outputStruct.Projection;
params.normDataCoeff=outputStruct.normDataCoeff;
params.substractMean=0;
params.sigma=0.02; %20 .020
params.delta=.050; %500
params.SDFWindowLength=0.020;%.020


coefficient = params.normDataCoeff;
auxCells =outputStruct.sd;
BinMode = outputStruct.BinMode;
ModalityFR = outputStruct.ModalityFr;
AvgCoef = outputStruct.AvgCoeff;
%-------------------------------------------%
%---------Get Behavior data-----------------%

if(strcmp(ModalityFR,'Average'))
    tempCells = reshape(auxCells,25,4);
    avgFr = cell(1);
    for i=1:4
      CondCells = tempCells(:,i)';
   avgFr{1,i}=mean(reshape(cell2mat(CondCells), [size(CondCells{1,1},1), size(CondCells{1,1},2), size(CondCells,2)]), 3);   

    end

    full_data = cat(2,avgFr{:});
    ToutputStruct.AvgSD = full_data ;
    
    if(params.normDataCoeff)
    normVector=max(full_data');
%     response = full_data';
%     normVector = (range(response) + 5);

    else
    normVector = ones(1,size(full_data,1));
    end
    %Coeff = 0; just transpose information;
    AvgfrData= bsxfun(@rdivide, full_data', normVector)*100;
    [Avgcoeff,scoresAvg, ~,~,explainedAvg]=pca(AvgfrData,'Economy',false);
    %------Average FirangRate Normalized - Proyections-----%
    avgFr = cellfun(@(x) bsxfun(@rdivide, x', normVector)*100 , avgFr, 'un' , 0); 
    avgFr = cellfun(@(x) fillmissing(x,'constant',0),avgFr,'un',0);
    AvgTraj = cellfun(@(x) getTrajectory(params,Avgcoeff,x'),avgFr,'un',0);
   
    Cells = cat(2,auxCells{:});
    ToutputStruct.AvgSD = Cells;
    clear dt normVector full_data
   
else
    Cells = cat(2,auxCells);
    AvgTraj = [];
    
end

%-----Cuttof Threshold - 5Hz-------%
if(outputStruct.threshold )
AVGFRTask = [avgFr{1,1}(11:98,:)',avgFr{1,2}(11:178,:)',avgFr{1,3}(11:98,:)',avgFr{1,4}(11:178,:)'];
% idx4Hz= find(any(AVGFRTask >= 10,2));
idx4Hz= find(mean(AVGFRTask,2) >= 6);

else
idx4Hz = [1:size(Cells,1)]';
end
%------------------------%


full_data = Cells(idx4Hz,:);
BinSize = 20;
%%Normalizing vector
if(params.normDataCoeff)
    normVector=max(full_data');
%     response = full_data';
%     normVector = (range(response) + 5);
%     respose_normalized = response/(range(response) + 5) %Soft normalization
      else
    normVector = ones(1,size(full_data,1));
end

a=1;

%Coeff = 1; Normalize
%Coeff = 0; just transpose information;
if(params.normDataCoeff)
    frData{a,1} = bsxfun(@rdivide, full_data', normVector);%*100;
else
    frData{a,1}=full_data';
end
frData{a,1}(isnan(frData{a,1}))=0;


CoeffPCA = outputStruct.coeffPCA;
dt = frData{1};

if(~isempty(CoeffPCA))
    coeff = outputStruct.coeffPCA;
else
%     [coeff,Sc, latent,~,explained]=pca(dt,'Economy',false);
    [coeff,Sc, latent,~,explained]=pca(dt);    

end


if(AvgCoef)
 Coeff =   Avgcoeff; 
%  [~,idx]= sort(Avgcoeff(:,1));
%  idx_300 = idx(1:500);
%  Coeff = Coeff(idx_300,idx_300);
%  normVector = normVector(idx_300);
%  idx4Hz = idx_300;
else
  Coeff =   coeff;
end


%Trials Projections
%Normalized FR
sd = cellfun(@(x) bsxfun(@rdivide, x(idx4Hz,:)', normVector)*100 , auxCells, 'un' , 0); 
sd = cellfun(@(x) fillmissing(x,'constant',0),sd,'un',0);

%Control - Task - Trajectories
traj = cellfun(@(x) getTrajectory(params,Coeff,x'),sd,'un',0); %coeff
%Control Trajectories
FRControl = outputStruct.FRControl;
if(~isempty(FRControl{1}))
FRControl = cellfun(@(x) bsxfun(@rdivide, x(idx4Hz,:)', normVector)*100 , FRControl, 'un' , 0); 
FRControl = cellfun(@(x) fillmissing(x,'constant',0),FRControl,'un',0);
traj30 = cellfun(@(x) getTrajectory(params,Coeff,x'),FRControl,'un',0); %coeff
else
traj30 = {[]};    
end






PrePostBin = 20;
SO = 4%5; %Serial Order

if(strcmp(outputStruct.BinMode,'Fixed')  )
   interval = [20 20 20 20]; 
else
interval = floor(outputStruct.interval / BinSize);
end
trials = outputStruct.trials;

ToutputStruct.expVar = explained;
ToutputStruct.expVarAvg = explainedAvg;

ToutputStruct.AvgTraj = AvgTraj;
ToutputStruct.traj = traj;
ToutputStruct.ControltrajFR = traj30;
ToutputStruct.Bin = (interval *4);
% ToutputStruct.Bin = [80 80 80 80];

ToutputStruct.BinMode = BinMode;
ToutputStruct.PrePostBin = PrePostBin;
ToutputStruct.Projection = params.projection;
ToutputStruct.trials = trials(1:length(interval));
ToutputStruct.SO = SO;%SO-1;
ToutputStruct.sd= sd;
ToutputStruct.Coeff= coeff;
ToutputStruct.AvgCoeff= Avgcoeff;
ToutputStruct.interval = outputStruct.interval;
ToutputStruct.Coefficient = coefficient;

outputStruct= ToutputStruct;












