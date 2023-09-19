function BinnedSDF = ComputeBinnesActiFromSDF(SSDF)

BinnSize = 25;%ms
TBins = numel(SSDF(1,:))/BinnSize;

tb = 1;
for bin = 1:TBins
    
    BinnedSDF(:,bin) = mean(SSDF(:,tb:tb+BinnSize-1)');    
    
    tb=tb+BinnSize;
end

% figure
%  imagesc(BinnedSDF);
 
%  figure
%  imagesc(SSDF);
 
 
 
 