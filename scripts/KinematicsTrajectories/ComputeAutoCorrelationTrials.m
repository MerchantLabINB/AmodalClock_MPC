function CorrePI = ComputeAutoCorrelationTrials(PI_meanTrials)


for k = 1:4
    for i=1:size(PI_meanTrials{1},1)
        PI_meanTrials{k}(isnan(PI_meanTrials{k})) = 0;
        [anglecorrPI(i,:),~] = autocorr(PI_meanTrials{k}(i,:));
    end
    CorrePI(k,:) = anglecorrPI(:,2);
end