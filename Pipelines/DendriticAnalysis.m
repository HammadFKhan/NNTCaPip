%% Dendrite Analysis
addpath(genpath('main'));
std_threshold = 6;
static_threshold = .01;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,length(Spikes));
cell_count = length(ROI);
time = time_adjust(size(DeltaFoverF,2),30.048);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,120);
end
factorCorrection = 5*floor(size(spikesTrig,2)/5); % Correct for frame size aquisition
[vectorized,sim_index] = cosine_similarity(spikesTrig(:,1:factorCorrection,1),5);
for i = 1:30
corr{i} = correlation_dice(spikesTrig(:,:,i));
Connected_ROI{i} = Connectivity_dice(corr{i}, ROI,0.3);
[NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
    ActivityCentroid,ActivityCentroidVariance]...
    = Network_Analysis(ROIcentroid,Connected_ROI{i});
end
%%
for i = 52
    thresh = std_threshold.*std(DeltaFoverF(i,:));
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