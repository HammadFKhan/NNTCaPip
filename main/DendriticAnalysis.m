%% Remove ROIs
if exist('badComponents','var') && ~exist('badComFlag','var')
    [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
        removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
    badComFlag = 1;
end
%%
ROIcentroid = [];ROInew = [];
for i = 1:length(ROI)
    trans = rand(1)*350;
    blah = vertcat(ROI{i}{:});
%     blah = ceil(blah+trans);
    polyin = polyshape(blah(:,1),blah(:,2));
    [x,y] = centroid(polyin);
    ROIcentroid(i,:) = [x y];
    for ii = 1:size(blah,1)
        ROInew{i,1}{ii,1} = blah(ii,:);
    end
end

%% Dendritic Coopertivity metric
% Generate a pseudo calcium fluorescence matrix (part of a covariance
% matrix) that generates traces during periods when a speific dendritic
% branch activates (can do 2-3 max) 
DenCoop = [];
count = 1;
for j = 1:length(IntraDendriticIndex) % intradendritic index format as each neuron is a cell with N array of dendrite pairs
    window = 300;
    idx = perms(IntraDendriticIndex{j});
    idx = unique(idx(:,[1 2]),'rows'); %take only the first two columns
    for i = 1:size(idx,1)
        thresh = 2.5.*std(DeltaFoverF(idx(i,1),:));
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
neuronID = 10;
primaryBranch = DenCoop(neuronID).branch{1};
sisterBranch = DenCoop(neuronID).branch{2};
figure,
for i = 1:size(primaryBranch,2)
subplot(ceil(sqrt(size(primaryBranch,2))),ceil(sqrt(size(primaryBranch,2))),i),...
    plot(primaryBranch(:,i)),axis off,ylim([0 6])
end

figure,
for i = 1:size(sisterBranch,2)
subplot(ceil(sqrt(size(sisterBranch,2))),ceil(sqrt(size(sisterBranch,2))),i)...
    ,plot(sisterBranch(:,i),'r'),axis off,ylim([0 6])
end

%% plot correlation of deltaF of branch pairs
primaryBranch = [];sisterBranch = [];
figure,
for i = 1:length(DenCoop)
    t1 = max(DenCoop(i).branch{1},[],1);
    t2 = max(DenCoop(i).branch{2},[],1);
    primaryBranch(i) = mean(t1);
    sisterBranch(i) = mean(t2);
    scatter(t1,t2,'fill','k'),hold on
end
axis([0 6 -0.1 6])
set(gca,'TickDir','out');
figure,boxplot(primaryBranch);box off,set(gca,'TickDir','out');ylim([0 5]);title('Primary Branch')
figure,boxplot(sisterBranch);box off,set(gca,'TickDir','out'),ylim([0 5]);title('Sister Branch')
% Coupling coefficient
figure,boxplot(abs(primaryBranch-sisterBranch));box off,set(gca,'TickDir','out'),ylim([0 5]);title('Coupling Coefficient')
%% Interdendritic activity (cross-population behavior)
addpath(genpath('main'));
std_threshold = 6;
static_threshold = .01;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
[coactive_cells,detected_spikes] = coactive_index(Spikes,floor(length(Spikes)*.05));

% Dendritic Population Coopertivity 
corr = correlation_dice(Spikes);
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);
corrEnsemble = correlation_dice(Ensemble.ensemble);
Connected_ROI = Connectivity_dice(corrEnsemble,0.15);
% Ensemble stats
Ensemble = ensembleMetric(Ensemble,AverageImage,ROIcentroid);
Ensemble = ensembleStat(Ensemble);
%% Plot
figure,imagesc(corr),colormap(jet),caxis([0 max(tril(corr,-1),[],'all')])
figure,Dendritic_Map(AverageImage,Connected_ROI,ROIcentroid,ROI,1,1),title('Pairwise Coopertivity')
figure,imagesc(corrEnsemble),colormap(jet),caxis([0 max(tril(corrEnsemble,-1),[],'all')]),colorbar
figure,Dendritic_Map(AverageImage,Connected_ROI,ROIcentroid,ROI,1,1)
