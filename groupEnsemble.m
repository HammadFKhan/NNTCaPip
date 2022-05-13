%% Single File Analysis
%% Set trial stats
ConensembleSize = [];
ConensembleRecruitment= [];
ConminRecruitment = [];
ConmaxRecruitment = [];
ConensembleNum = [];
ConeCosine = [];
% W6simIndex = [];
ConSDactivity = [];
% load batch directory
pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for idx = 1:L
load([pathname directory(idx).name])
addpath(genpath('main'));
std_threshold = 4;
static_threshold = .2;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,length(Spikes));
cell_count = length(ROI);
time = time_adjust(num_images,15);
% calcium_avg = STA(DeltaFoverF,Spikes,std_threshold,5);
bin = 20; %Vector sizes for similarity indexing (Num frames should be devisable by this)           

% Spikes_shuffled = tempShuffle(Spikes,10000);
% Event_shuffled = spatialShuffle(Spikes,10000);
% surrogate = 10;
% Total_shuffled = allShuffle(Spikes,10000);
% [shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
% shuff_corr = correlation_dice(Event_shuffled);
% [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
% shufsim_index = shufsim_index-mean(mean(shufsim_index,2));
% factorCorrection = 10*floor(size(Spikes,2)/10); % Correct for frame size aquisition
% [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),10);
corr = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr, ROI);
[NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
    ActivityCentroid,ActivityCentroidVariance]...
    = Network_Analysis(ROIcentroid,Connected_ROI);
% Ensemble Analysis
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROI,ROIcentroid);
% Ensemble = ensembleNetworks(Ensemble);

% %% Plot Ensemble
% ensembleVid(Ensemble,AverageImage,ROIcentroid,files);
% 
% % Combine Maps
% figure,imagesc(interp2(Ensemble.sim_index,2)),colormap(jet),caxis([0.13 0.4])
% K = (1/25)*ones(5);
% figure,imagesc(interp2(conv2(Ensemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .15])

%
    % Ensemble stats for trial processing
    [ensembleSize,ensembleNum,ensembleRecruitment,minRecruitment,maxRecruitment] = ensembleStat(Ensemble,ROI);
    ConeCosine = vertcat(ConeCosine,mean(Ensemble.sim_index,2));
    ConensembleSize = vertcat(ConensembleSize,ensembleSize');
    ConensembleNum = vertcat(ConensembleNum,ensembleNum);
    ConensembleRecruitment= vertcat(ConensembleRecruitment,ensembleRecruitment);
    ConminRecruitment = vertcat(ConminRecruitment,minRecruitment);
    ConmaxRecruitment = vertcat(ConmaxRecruitment,maxRecruitment);
%     meanSim = mean(sim_index);
%     W6simIndex = vertcat(W6simIndex,meanSim);
    meanSDactivity = mean(corr,2);
    ConSDactivity = vertcat(ConSDactivity,meanSDactivity);
    %%
    clearvars -except ConensembleSize ConensembleRecruitment ConminRecruitment ConmaxRecruitment ConensembleNum...
        ConSDactivity ConeCosine ConsimIndex ConSDactivity files pathname directory L
%         W8ensembleSize W8ensembleRecruitment W8minRecruitment W8maxRecruitment...
%         ConensembleSize ConensembleRecruitment ConminRecruitment ConmaxRecruitment
end