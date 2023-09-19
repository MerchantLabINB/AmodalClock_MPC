function   [AngleTaps,DistProyeArr,MargangleArr,VarPosition,StartingPointArr] = KinematicsParamsProyectionTrials(LinePoint,PointMod,FitPoints,scores,scoresM,PlanePoints,binarrTaps,Amplitude)


nCond = 4;
numPC = 3;
numPcVar = 1;
for numPoint = LinePoint%1: size(FitPoints,1)% 224 358 116
    TapsAvgPoint = FitPoints(numPoint,1:3);
    %Distance Implmentation

    for i = 1:nCond%nCond:-1:1
        DistProye = [];
        Margangle_aux = [];
        StartingPoint = [];
        VaraibilityfromPos = [];
        for numtrials = 1:25
            for SO = 1:4%4
                ntt = size(scores{numtrials,i},1);
                for j =1:ntt
                    %% Match time points by distance (Speed, Distance)
                    ep_data = scores{numtrials,i}(j,1:3,SO);
                    %                     ep_dataM = scoresM{numtrials,i}(j,1,SO);
                    ep_dataM = scoresM{numtrials,i}(j,1:3,SO);

                    ModPoint = PlanePoints(PointMod,:);
                    distances = sqrt(sum((TapsAvgPoint' - ...
                        scores{numtrials,i}(j,1:3,SO)').^2));

                    DistProye( numtrials,binarrTaps{1,i}(j,SO))= distances;
                    StartingPoint(numtrials,binarrTaps{1,i}(j,SO),:) = ep_dataM-ModPoint;%ep_data-ModPoint;%sum(ep_data-TapsAvgPoint);
                    aux = norm(TapsAvgPoint) * norm(ep_data);
                    Margangle( binarrTaps{1,i}(j,SO),i) = acos((dot(TapsAvgPoint,ep_data) /aux ));
                    Margangle_aux(numtrials, binarrTaps{1,i}(j,SO)) = acos((dot(TapsAvgPoint,ep_data) /aux ));
                    %% Variability of trajectories
                end
                [val,idxTaps]= min(Margangle(binarrTaps{1,i}(:,SO),i));

                if(idxTaps == 1)
                    [val,indexes] = sort(Margangle(binarrTaps{1,i}(:,SO),i));
                    idxTaps = indexes(2);
                    val = val(2);
                end
                addVal = binarrTaps{1,i}(1,SO) -1;
                idxTap(numtrials,SO,i) = idxTaps + addVal ;
                aux_Taps(numtrials,SO) = (idxTaps + addVal)*20;
                VaraibilityfromPos(numtrials,SO) = std(StartingPoint(numtrials,binarrTaps{1,i}(:,SO),numPcVar));
                %                 VaraibilityfromPos(numtrials,SO) = std(StartingPoint(numtrials,binarrTaps{1,i}(:,SO)));
                minAngle(numtrials,SO) = val;
            end
        end

        AngleTaps{i} = [ones(25,1),aux_Taps];
                DistProyeArr{i} = DistProye- DistProye(:,1);
%         DistProyeArr{i} = DistProye;


        StartingPointArr{i} = squeeze(StartingPoint(:,:,numPcVar));
        PhaseAngle{i} = minAngle;
        VarPosition{i} = VaraibilityfromPos;
        MargangleArr{i} = Margangle_aux;
    end

end
