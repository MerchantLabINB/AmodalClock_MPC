function [meanAutoCorrAll,meanAutoCorrFirstLag,autocorrAllTrials] = lagAutoCorrelation(int)
for idxInt = 1:size(int,2)
    for idxTrial = 1:size(int{idxInt},1)
        
        autocorrAllTrials{1,idxInt}(idxTrial,:) = autocorr(int{idxInt}(idxTrial,:));
        
    end
    meanAutoCorrAll{1,idxInt} = mean(autocorrAllTrials{1,idxInt},1);
     meanAutoCorrFirstLag{1,idxInt} = meanAutoCorrAll{1,idxInt}(1,2);
end