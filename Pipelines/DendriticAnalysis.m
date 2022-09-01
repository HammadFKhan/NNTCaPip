%% Dendrite Analysisaddpath(genpath('main'));
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
figure,
for i = 1:30
%     title(['Trial ' num2str(i)])
    subplot(6,5,i),imagesc(sim_index(:,:,i)),colormap(hot),axis off;
end

figure,
for i = 1:30
    subplot(6,5,i),imagesc(corr{i});colormap(jet),axis off
end
figure,
for i = 1:30
subplot(6,5,i),Cell_Map_Dice(AverageImage,Connected_ROI{i},ROIcentroid,NodeSize,EdgeSize)
end
%%
addpath('Figures');
figure('Name','DeltaF/F'); stack_plot(DeltaFoverF,1.5,15); 
figure('Name','Convolved Spikes'); plot(dDeltaFoverF');
figure('Name','Threshold Detection');DeltaFoverFplotter(dDeltaFoverF,std_threshold,static_threshold)
figure('Name','Spike Plot'); Show_Spikes(Spikes);
% figure('Name','Temporal Shuffled Spike Plot'); shuffledTspikePlot = Show_Spikes(Spikes_shuffled);
% figure('Name','Event Shuffled Spike Plot'); shuffledEspikePlot = Show_Spikes(Event_shuffled);
% figure('Name','Total Shuffled Spike Plot'); shuffledAspikePlot = Show_Spikes(Total_shuffled);
figure('Name','Fluorescence Map'); spike_map(DeltaFoverF);caxis([0 1]),set(gcf,'Position',[100 100 400 400])
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Coactivity Index'); B = bar(coactive_cells,4);ax = gca;ax.TickDir = 'out';ax.Box = 'off';
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.2]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Dice-Similarity Index');h = htmp(shuff_corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([.7 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Cosine-Similarity Index'); h = htmp(shufsim_index);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Network Map'); NodeSize = 2;EdgeSize = 2;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize)
