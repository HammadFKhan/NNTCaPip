%% Run the CaimAn program
% Github Version 4.1

%####Batch Analysis for ROI Extaction####
addpath(genpath('main'));
clear
clc
close all;
tic
CaImAn_Conversion
disp('Batch Analysis complete')
toc
%% Batch Data Analysis
addpath(genpath('main'));
batchData = batchData_Analysis;
%% CaimAn Single File ROI Extraction 
clear
clc
close all;
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
tic
Start_CaImAn
toc
%% CaimAn Single File ROI Extraction from Dendrites
clear
clc
close all;
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
tic
Start_CaImAn_Dendrites
toc
%% Single File Analysis
set(0,'DefaultFigureWindowStyle','docked')
addpath(genpath('main'));
cell_count = length(ROI);
time = time_adjust(num_images,15);
std_threshold = 5;
static_threshold = .2;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
corr = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr, ROI);
corr_jaccard = correlation_jaccard(Spikes);
[coactive_cells,detected_spikes] = coactive_index(Spikes);
calcium_avg = STA(DeltaFoverF,Spikes,std_threshold,5);
bin = 20; %Vector sizes for similarity indexing (Num frames should be devisable by this)           

Spikes_shuffled = tempShuffle(Spikes,num_images,cell_count);
Event_shuffled = spatialShuffle(Spikes,num_images,cell_count);
surrogate = 10;
Total_shuffled = allShuffle(Spikes,num_images,cell_count,surrogate);
shuff_corr = correlation_dice(Total_shuffled);

[vectorized,sim_index] = cosine_similarity(Spikes,bin);
[shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
shufsim_index = shufsim_index-mean(mean(shufsim_index,2));
[NumActiveNodes,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
    ActivityCentroid,ActivityCentroidVariance]...
    = Network_Analysis(ROIcentroid,Connected_ROI);

a = mean(mean(sim_index,2));
b = mean(abs(mean(shufsim_index,2)));
%% Trial by Trial analysis
addpath(genpath('Figures'));
[batchSpikes,batch_corr] = TrialByTrial(batchData([1,4,5,6,7,8])); % Function call
bin = 20;
[vectorized,sim_index] = cosine_similarity(batchSpikes,bin);
[z,mu,sigma] = zscore(sim_index);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);
% caxis([.5 0.9]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc', 'Figures/CSTrials.eps', '-r250');
 figure('Name','Dice Correlation')
 for i = 1:size(batch_corr,3)
     subplot(2,3,i),h = htmp(batch_corr(:,:,i),20);caxis([0 0.4]);
 end
%% Plot all the Figures
addpath('Figures');
figure('Name','DeltaF/F'); stack_plot(DeltaFoverF,0.8); 
% figure('Name','Convolved Spikes'); plot(dDeltaFoverF');
% figure('Name','Threshold Detection');DeltaFoverFplotter(dDeltaFoverF,std_threshold,static_threshold)
figure('Name','Spike Plot'); spikePlot = Show_Spikes(Spikes);
figure('Name','Temporal Shuffled Spike Plot'); shuffledTspikePlot = Show_Spikes(Spikes_shuffled);
figure('Name','Event Shuffled Spike Plot'); shuffledEspikePlot = Show_Spikes(Event_shuffled);
figure('Name','Total Shuffled Spike Plot'); shuffledAspikePlot = Show_Spikes(Total_shuffled);
figure('Name','Fluorescence Map'); spikeImage = spike_map(DeltaFoverF,time);caxis([0 .4]);...
    print(gcf,'-painters','-depsc', 'Figures/Fluormap.eps', '-r250');
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
    print(gcf,'-painters','-depsc', 'Figures/Population.eps', '-r250');
figure('Name','Coactivity Index'); B = bar(time,coactive_cells);...
    ax = gca;ax.TickDir = 'out';ax.Box = 'off';axis off;...
    print(gcf,'-painters','-depsc', 'Figures/Coactive.eps', '-r250');
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc', 'Figures/DSC.eps', '-r250');
figure('Name','Shuffled Dice-Similarity Index');h = htmp(shuff_corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc', 'Figures/sDCS.eps', '-r250');
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.35 .9]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc', 'Figures/CS.eps', '-r250');
figure('Name','Shuffled Cosine-Similarity Index'); h = htmp(shufsim_index);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc', 'Figures/sCS.eps', '-r250');
figure('Name','Network Map'); NodeSize = 6;EdgeSize = 2;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize)

%% Rotary Encoder
figure('Name','Pulse Data');plot(encoder_data.rotate_pulse);
figure('Name','Angular Distance');bar(encoder_data.ang_distance);
figure('Name','Angular Velocity');bar(encoder_data.ang_velocity,'FaceColor',[.16 .835 .384],'EdgeColor','none');
figure('Name','Avg. Angular Velocity');avgV = movmean(encoder_data.ang_velocity,2);bar(avgV,'FaceColor',[.16 .835 .384],'EdgeColor','none');

%%
image_movie = mat2gray(Image_Stack);
implay(image_movie);

% Have you tried using Multidimensional Scaling (MDS) to emebed the
% centroids in a 2 dimensional space for visualization?

% This should visualize how the centroids related to each other. You could�
% also then compute the Delauney Triangulation of the projected graph, to
% identify neighbors.
