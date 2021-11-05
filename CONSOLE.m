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
std_threshold = 2.5;
static_threshold = .2;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
[coactive_cells,detected_spikes] = coactive_index(Spikes,length(Spikes));
calcium_avg = STA(DeltaFoverF,Spikes,std_threshold,5);
bin = 10; %Vector sizes for similarity indexing (Num frames should be devisable by this)           

Spikes_shuffled = tempShuffle(Spikes,10000);
% Event_shuffled = spatialShuffle(Spikes,10000);
% surrogate = 10;
% Total_shuffled = allShuffle(Spikes,10000);
[shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
shuff_corr = correlation_dice(Event_shuffled);
% [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
% shufsim_index = shufsim_index-mean(mean(shufsim_index,2));

%%
[vectorized,sim_index] = cosine_similarity(Spikes(:,1:4180),bin);
corr = correlation_dice(Spikes);
Connected_ROI = Connectivity_dice(corr, ROI);
[NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
    ActivityCentroid,ActivityCentroidVariance]...
    = Network_Analysis(ROIcentroid,Connected_ROI);
%% Rank Correlation
rank = mean(corr);
[~,idx_sort] = sort(rank);
rankCorr = corr(idx_sort,:);
%%
idx = find(trig==1); %Find when the trigger is high
onset = find(diff(idx)>1)+1; %Find when the trigger ends
events(:,1) = idx(:,1)/20000; %Take only the first value of the trigger event
for i = 1:length(onset)
    events(:,i+1) = idx(:,onset(:,i))/20000; % Map the onset of the pulse in seconds
end
%%
Peristim{1} = Spikes(:,1:120);
for i = 1:floor(4182/120)-1 % Extract calcium data +-2s to stim (120 frames per trial) number should equal trials
    Peristim{i+1} = Spikes(:,(i*120)+1:((i+1)*120));
end
Peristim{ceil(4182/120)} = Spikes(:,(i+1)*120:end); %Account for differences at the end of the trail
%%
for i = 1:length(ensembleCan)
    ensemble(:,i) = stimTrials(:,ensembleCan(i));
end

%% Ensemble Analysis
% win = find(t1==1);
% UDS = find(diff(win)>1)
% fh1 = figure;
% fh2 = figure;
countS = 1;
j = 6;
bins = discretize(1:length(Peristim{j}(1,:)),60); %~13 frames
for j = 3
    for i = 1:8
        disp(['Detecting Ensembles in bin: '  num2str(i)])
        bins = discretize(1:length(Peristim{j}(1,:)),8); %~13 frames
        corr = correlation_dice(Peristim{j});
        %             avgCorr(:,i) = mean((mean(corr)));
        %             [coactive_cells,detected_spikes] = coactive_index(lateSpike{j},30);
                        Connected_ROI{i} = Connectivity_dice(corr, ROI);
                        [NumActiveNodes(j,i),NodeList,NumNodes(j,i),NumEdges(j,i),SpatialCentroid,SpatialCentroidVariance,...
                            ActivityCentroid,ActivityCentroidVariance]...
                            = Network_Analysis(ROIcentroid,Connected_ROI);
        %         Extract active nodes to graph
        %                 id = find(Connected_ROI(:,3)>0.6);
        %                 ensembleNodes = [Connected_ROI(id,1) Connected_ROI(id,2)];
        %                 strongConnections(j,i) = length(find(Connected_ROI(:,3)>0.5));
        %                 weakConnections(j,i) = length(find(Connected_ROI(:,3)<0.5));
        %     figure(1),h = htmp(corr);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
        %                     figure(fh1), subplot(2,4,i),NodeSize = 6;EdgeSize = 3;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize);
        %             title(['Time: ' num2str(j) ' seconds']);
        %                     count = count+1;
        %         corr = correlation_dice(ensemble(:,80:100));
        %             Connected_ROI = Connectivity_dice(corr, ROI);
        %             [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
        %                 ActivityCentroid,ActivityCentroidVariance]...
        %                 = Network_Analysis(ROIcentroid,Connected_ROI);
        %             hold on;
        %             figure(2),
        %             EnsembleMap(AverageImage,ROIcentroid,NodeList,NodeSize)
        %     drawnow
        %             hold off;
    end
end
%% %% Ensemble Analysis
% win = find(t1==1);
% UDS = find(diff(win)>1)
% fh1 = figure;
% fh2 = figure;
countS = 1;
j = 6;
bins = discretize(1:length(Peristim{j}(1,:)),60); %~13 frames
for j = 1:16
    disp(['Detecting Ensembles in bin: '  num2str(i)])
    bins = discretize(1:length(weaklateSpike{j}(1,:)),8); %~13 frames
    corr = correlation_dice(weaklateSpike{j});
    %             avgCorr(:,i) = mean((mean(corr)));
    %             [coactive_cells,detected_spikes] = coactive_index(lateSpike{j},30);
    %                 Connected_ROI{i} = Connectivity_dice(corr, ROI);
    %                 [NumActiveNodes(j,i),NodeList,NumNodes(j,i),NumEdges(j,i),SpatialCentroid,SpatialCentroidVariance,...
    %                     ActivityCentroid,ActivityCentroidVariance]...
    %                     = Network_Analysis(ROIcentroid,Connected_ROI);
    %         Extract active nodes to graph
    %                 id = find(Connected_ROI(:,3)>0.6);
    %                 ensembleNodes = [Connected_ROI(id,1) Connected_ROI(id,2)];
    %                 strongConnections(j,i) = length(find(Connected_ROI(:,3)>0.5));
    %                 weakConnections(j,i) = length(find(Connected_ROI(:,3)<0.5));
%     figure(1),h = htmp(corr);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
    %                     figure(fh1), subplot(2,4,i),NodeSize = 6;EdgeSize = 3;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize);
%             title(['Time: ' num2str(j) ' seconds']);
    %                     count = count+1;
    %         corr = correlation_dice(ensemble(:,80:100));
            Connected_ROI = Connectivity_dice(corr, ROI);
            [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
                ActivityCentroid,ActivityCentroidVariance]...
                = Network_Analysis(ROIcentroid,Connected_ROI);
            hold on;
            figure(2),
            EnsembleMap(AverageImage,ROIcentroid,NodeList,NodeSize)
    drawnow
            hold off;
end
%%
for j = 1:length(Peristim)
    for i = 2:4
        stack{j,i-1} = Peristim{j}(:,find(bins==i));
    end
    responseCorrelation{j} = horzcat(stack{j,:});
end
responseCorrelation = cat(2,responseCorrelation{:});
%%
corr = correlation_dice(ensemble(:,80:100));
Connected_ROI = Connectivity_dice(corr, ROI);
[NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
                    ActivityCentroid,ActivityCentroidVariance]...
                    = Network_Analysis(ROIcentroid,Connected_ROI);
figure,
EnsembleMap(AverageImage,ROIcentroid,NodeList,NodeSize)
%% Concatenate no stim trails
for i = 1:length(allTrials)
    nostimTrials{i} = Peristim{allTrials(i,1)}(:,end-44:end);
end
nostimTrials = horzcat(nostimTrials{:});
%% Population Activity before and after stim
[vectorized,sim_index] = cosine_similarity(Spikes(:,1:4175),25);
% figure,imagesc(vectorized),colormap(jet),colorbar
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.35 .9]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
[vectorized,sim_index] = cosine_similarity(nostimTrials,25);
figure,imagesc(vectorized),colormap(jet),colorbar
% figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.35 .9]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
% figure('Name','Fluorescence Map'); spikeImage = spike_map(stimTrials,time);caxis([0 1]);
% figure('Name','Fluorescence Map'); spikeImage = spike_map(nostimTrials,time);caxis([0 .1]);
%%
corr = correlation_dice(stimTrials);
Connected_ROI = Connectivity_dice(corr, ROI);
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.3]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
[NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
    ActivityCentroid,ActivityCentroidVariance]...
    = Network_Analysis(ROIcentroid,Connected_ROI);
figure,
EnsembleMap(AverageImage,ROIcentroid,NodeList,5);
%%
ensembleCan = find(coactive_cells>0.13);
for i = 1:length(ensembleCan)
ensemble(:,i) = Spikes(:,ensembleCan(i));
end
%% 
[vectorized,sim_index] = cosine_similarity(ensemble(:,1:190),2);
svd_analysis(sim_index)

%%
% Strong vs Weak input response
% Checks the stong vs weak connection 1-2 seconds after pole swing
countS = 1;
countW = 1;
for i = 1:length(LateSpikeTrials)
        weakLateResponse(i,:) = [strongConnections(CountS,:), weakConnections(WeakLateSpikeTrials(i,:),:)];
end

%% 
trial = 1:12;
figure
plot(mean(LateResponse(trial,1:8),1)), hold on;
plot(mean(LateResponse(trial,9:16,1)))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
trial = 1:16;
figure
plot(mean(weakLateResponse(trial,1:8),1)), hold on;
plot(mean(weakLateResponse(trial,9:16,1)))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
%%
figure,plot(strongConnections(:,8))
hold on,plot(weakConnections(:,8))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);

%%
figure
for i = 1:size(strongConnections,1)
    plot(strongConnections(i,:)','Color',[0.2 0.2 0.2 i/25]), hold on;
end
figure
for i = 1:size(strongConnections,1)
    plot(weakConnections(i,:)','Color',[0.2 0.2 0.2 i/25]), hold on;
end
ylim([0 200])
%% Trial by Trial analysis ##Only use with batch processed files##
addpath(genpath('Figures'));
[batchSpikes,batch_corr] = TrialByTrial(batchData([1,2,4])); % Function call
bin = 20;
[vectorized,sim_index] = cosine_similarity(batchSpikes,bin);
[z,mu,sigma] = zscore(sim_index);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index,100);
caxis([0 0.7]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);...
     print(gcf,'-painters','-depsc2', 'Figures/CSTrials.eps', '-r250');
 figure('Name','Dice Correlation')
 for i = 1:size(batch_corr,3)
     subplot(2,3,i),h = htmp(batch_corr(:,:,i),20);caxis([0 0.4]);
 end
%% Plot all the Figures
addpath('Figures');
figure('Name','DeltaF/F'); stack_plot(DeltaFoverF,0.8,2); 
figure('Name','Convolved Spikes'); plot(dDeltaFoverF');
figure('Name','Threshold Detection');DeltaFoverFplotter(dDeltaFoverF,std_threshold,static_threshold)
figure('Name','Spike Plot'); spikePlot = Show_Spikes(Spikes);
% figure('Name','Temporal Shuffled Spike Plot'); shuffledTspikePlot = Show_Spikes(Spikes_shuffled);
% figure('Name','Event Shuffled Spike Plot'); shuffledEspikePlot = Show_Spikes(Event_shuffled);
% figure('Name','Total Shuffled Spike Plot'); shuffledAspikePlot = Show_Spikes(Total_shuffled);
figure('Name','Fluorescence Map'); spikeImage = spike_map(DeltaFoverF,time);caxis([0 .45]);
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Coactivity Index'); B = bar(coactive_cells);ax = gca;ax.TickDir = 'out';ax.Box = 'off';
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.3]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Dice-Similarity Index');h = htmp(shuff_corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.35 .9]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Cosine-Similarity Index'); h = htmp(shufsim_index);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Network Map'); NodeSize = 2;EdgeSize = 2;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize)

%% Rotary Encoder
figure('Name','Pulse Data');plot(encoder_data.rotate_pulse);
figure('Name','Angular Distance');bar(encoder_data.ang_distance);
figure('Name','Angular Velocity');bar(encoder_data.ang_velocity,'FaceColor',[.16 .835 .384],'EdgeColor','none');
figure('Name','Avg. Angular Velocity');avgV = movmean(encoder_data.ang_velocity,2);bar(avgV,'FaceColor',[.16 .835 .384],'EdgeColor','none');

%%
image_movie = mat2gray(Image_Stack);
implay(image_movie);
%%
plot_raster(1:120,Spikes(5,1:120))
% Have you tried using Multidimensional Scaling (MDS) to emebed the
% centroids in a 2 dimensional space for visualization?

% This should visualize how the centroids related to each other. You could�
% also then compute the Delauney Triangulation of the projected graph, to
% identify neighbors.
function [tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(featureData,navgtarget)
featureData = double(featureData');
covmatrix = (featureData'*featureData);
covmatrix = covmatrix/size(featureData,1);
figure();
imagesc(covmatrix);
colormap(jet);
colorbar;
[V,D] = eig(covmatrix);
q(:,1) = V(:,size(featureData,2));
q(:,2) = V(:,size(featureData,2)-1);
q(:,3) = V(:,size(featureData,2)-3);
figure();
plot(q);
ylabel('Voltage (\mu V)')
xlabel('Time');

tProjq1 = featureData(1:navgtarget,:)*q(:,1);
tProjq2 = featureData(1:navgtarget,:)*q(:,2);
tProjq3 = featureData(1:navgtarget,:)*q(:,3);
figure()
scatter3(tProjq1,tProjq2,tProjq3,200,'b.'); hold on;
uProjq1 = featureData(navgtarget+1:end,:)*q(:,1);
uProjq2 = featureData(navgtarget+1:end,:)*q(:,2);
uProjq3 = featureData(navgtarget+1:end,:)*q(:,3);
scatter3(uProjq1,uProjq2,uProjq3,200,'r.');
ylabel('bk')
xlabel('ak');
end
