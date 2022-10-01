%% Dendritic Coopertivity metric
% Generate a pseudo calcium fluorescence matrix (part of a covariance
% matrix) that generates traces during periods when a speific dendritic
% branch activates (can do 2-3 max) 
ROIcentroid = [];
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    polyin = polyshape(blah(:,1),blah(:,2));
    [x,y] = centroid(polyin);
    ROIcentroid(i,:) = [x y];
end
%%
count = 1;
for j = 1:length(IntraDendriticIndex)
    window = 300;
    idx = perms(IntraDendriticIndex{j});
    idx = unique(idx(:,[1 2]),'rows'); %take only the first two columns
    for i = 1:size(idx,1)
        thresh = 2.*std(DeltaFoverF(idx(i,1),:));
        [pks,locs] = findpeaks(DeltaFoverF(idx(i,1),:),'MinPeakHeight',thresh);
        spikeCount = 0;
        if(~isempty(pks))
            for ii = 1:length(locs)
                thisLoc = locs(ii);
                if(thisLoc-window>1 && thisLoc+window<length(DeltaFoverF))%%peak can't occur at the very end or very beginning of the data set
                    spikeCount = spikeCount+1;
                    %extract a 4ms window around the spike peak
                    DenCoop(count).branch{1}(:,spikeCount) = DeltaFoverF(idx(i,1),thisLoc-(window/2):thisLoc+window); % take the primary branch
                    DenCoop(count).branch{2}(:,spikeCount) = DeltaFoverF(idx(i,2),thisLoc-(window/2):thisLoc+window); % compare to sister branch
                    DenCoop(count).branchLabel = idx(i,:);
                end
            end
        end
        count = count+1;
    end
end
%%
figure,
for i = 1:15
subplot(3,5,i),plot(DenCoop(13).branch{2}(:,i)),axis off,ylim([0 6])
end
%% plot correlation of deltaF of branch pairs
figure,
for i = 1:length(DenCoop)
    t1 = max(DenCoop(i).branch{1},[],1);
    primaryBranch(i) = mean(t1);
    sisterBranch(i) = mean(t2);
    t2 = max(DenCoop(i).branch{2},[],1);
    scatter(t1,t2,'fill','k'),hold on
end
axis([0 6 -0.1 6])
set(gca,'TickDir','out');
%% Interdendritic activity (cross-population behavior)
factorCorrection = 10*floor(size(Spikes,2)/10); % Correct for frame size aquisition
% [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),10);
corr = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr, ROI,0.15);
% [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
%     ActivityCentroid,ActivityCentroidVariance]...
%     = Network_Analysis(ROIcentroid,Connected_ROI);
figure,Dendritic_Map(AverageImage,Connected_ROI,ROIcentroid,ROI,0,1)