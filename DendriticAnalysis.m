%% Dendritic Coopertivity metric
% Generate a pseudo calcium fluorescence matrix (part of a covariance
% matrix) that generates traces during periods when a speific dendritic
% branch activates (can do 2-3 max) 
idx = 1;
for i = [105 106]
    thresh = 2.*std(DeltaFoverF(i,:));
    [pks,locs] = findpeaks(DeltaFoverF(i,:),'MinPeakHeight',thresh);
    spikeCount = 0;
    if(~isempty(pks))
        for ii = 1:length(locs)
            thisLoc = locs(ii);
            if(thisLoc-window>1 && thisLoc+window<length(DeltaFoverF))%%peak can't occur at the very end or very beginning of the data set
                spikeCount = spikeCount+1;
                %extract a 4ms window around the spike peak
                 ROI2{2}(:,spikeCount) = DeltaFoverF(105,thisLoc-(window/2):thisLoc+window);
            end
        end
    end
    idx = idx+1;
end

%%
figure,
for i = 1:18
subplot(5,4,i),plot(ROI1{2}(:,i)),axis off,ylim([0 3.5])
end