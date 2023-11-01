function calcium_avg = STA(DeltaFoverF,std_threshold,window)

%% find all spikes and align to the peak
%find negative peaks in filtered data with an amplitude of at least 4 std and 1ms  in
%separation
calcium_avg = [];
%for each detected spike, extract a window 2ms before and after the peak
for i = 1:length(DeltaFoverF(:,1))
    thresh = std_threshold*std(DeltaFoverF(i,:));
    [pks,locs] = findpeaks(DeltaFoverF(i,:),'MinPeakHeight',thresh);
    spikeCount = 0;
    if(~isempty(pks))
        for ii = 1:length(locs)
            thisLoc = locs(ii);
            if(thisLoc-window>1 && thisLoc+window<length(DeltaFoverF))%%peak can't occur at the very end or very beginning of the data set
                spikeCount = spikeCount+1;
                %extract a 4ms window around the spike peak
                 calcium_avg(:,spikeCount) = DeltaFoverF(thisLoc-(window/2):thisLoc+window);
            end
        end
    end
end

end