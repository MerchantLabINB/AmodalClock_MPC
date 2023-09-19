function   [DistProye,Margangle,StartingPoint,idxTap,VarPosition] = KinematicsParamsProyection(LinePoint,PointMod,FitPoints,scores,scoresM,PlanePoints,binarrTaps)
nCond = 4;
numPcVar = 3;

for numPoint = LinePoint
    TapsAvgPoint = FitPoints(numPoint,1:3);
    for SO = 1:4

        %Distance Implmentation
        %Avg Taps
        for i = 1:nCond
            ntt = size(scores{1,i},1);
            for j =1:ntt

                %% Match time points by distance (Speed, Distance)
                ep_data = scores{1,i}(j,1:3,SO);
                ep_dataM = scoresM{1,i}(j,1:3,SO);
                ModPoint = PlanePoints(PointMod,:);

                ep_data_marg = ep_data./ norm(ep_data);
                distances = sqrt(sum((TapsAvgPoint' - ...
                    scores{1,i}(j,1:3,SO)').^2));

                minDistTap( binarrTaps{1,i}(j,SO),i)= distances;
                DistProye( binarrTaps{1,i}(j,SO),i)= distances;
                StartingPoint(binarrTaps{1,i}(j,SO),i,:) = ModPoint- ep_dataM;%ep_data-ModPoint;%sum(ep_data-TapsAvgPoint);
                aux = norm(TapsAvgPoint) * norm(ep_data);
                Margangle( binarrTaps{1,i}(j,SO),i) = acos((dot(TapsAvgPoint,ep_data) /aux ));


            end

            [val,idxTaps]= min(Margangle(binarrTaps{1,i}(1:end,SO),i));

            if(idxTaps < 3)
                a = 1;
                iter = 0;
                while(idxTaps <= 3)
                    [val,idxTaps]= min(Margangle(binarrTaps{1,i}(a:end,SO),i));
                    a = a+1;
                    iter = iter+1;
                end
                idxTaps = idxTaps + (iter-1);
            end

            addVal = binarrTaps{1,i}(1,SO) -1;
            idxTap(i,SO) = idxTaps + addVal ;
            VaraibilityfromPos(i,SO) = std(StartingPoint(binarrTaps{1,i}(:,SO),i,numPcVar));
        end

    end
    VarPosition  = num2cell(VaraibilityfromPos,2)';
end
idxTap = [ones(4,1),idxTap *20]; %miliseconds