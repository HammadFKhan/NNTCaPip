%% NNT Calcium Pipeline
% Github Version 4.1

% CaimAn Single File ROI Extraction 
clear
clc
close all; 
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
tic
Start_CaImAn
toc
%% CaimAn Single File ROI Extraction (Large Dataset)
clear
clc
close all;
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
global nam
global memfig
batch = 0;
if batch == 1
    pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
    savepathname = strcat(uigetdir(pwd,'Output Directory'),'\');
    directory = dir(pathname);
    L = length(directory);
    for idx = 3:L
        clear files AverageImage num_images...
            DeltaFoverF dDeltaFoverF ROIcentroid ROI Noise_Power
        filename = directory(idx).name
        nam = strcat(pathname,filename);
        tic
        Start_MemMap_CaImAn
        toc
        savepath = strcat(savepathname,filename,'.mat');
        save(savepath,'files','AverageImage','num_images',...
            'DeltaFoverF','dDeltaFoverF','ROIcentroid','ROI','Noise_Power');
        try
            savepathfig = strcat(savepathname,filename(1:end-4),'.png');
            saveas(memfig,savepathfig);
        catch ME
            warning('Contour figure not saved')
            continue
        end
        disp('Saved!')
    end
else
    nam = '';
    tic
    Start_MemMap_CaImAn
    toc
end
%% Analysis
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
time = time_adjust(size(DeltaFoverF,2),30);
calcium_avg = STA(DeltaFoverF,Spikes,std_threshold,5);

% Perform shuffling and pairwise if data is small enough
if num_images<2000    
%     Spikes_shuffled = tempShuffle(Spikes,1000);
%     Event_shuffled = spatialShuffle(Spikes,1000);
%     surrogate = 10;
%     Total_shuffled = allShuffle(Spikes,1000);
%     [shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
%     shuff_corr = correlation_dice(Event_shuffled);
%     [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
%     shufsim_index = shufsim_index-mean(mean(shufsim_index,2));
    factorCorrection = 10*floor(size(Spikes,2)/10); % Correct for frame size aquisition
    [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),10);
    corr = correlation_dice(Spikes);
    Connected_ROI = Connectivity_dice(corr, ROI,0.3);
    [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
        ActivityCentroid,ActivityCentroidVariance]...
        = Network_Analysis(ROIcentroid,Connected_ROI);
end

%% Pairwise synchronization as velocity
Velocity = encoderVelocity(VR_data); % Bin width of 4 (ie. average every 4 readings)
% Extract window during velocity threshold
thresh = find(abs(Velocity(:,2))>1);
thresh = thresh*4; %multiply by bin width of the encoder time
% Now match threshold to the calcium time by window
for i = 1:length(thresh)
    CaTime = find(round(time)==thresh(i));
    syncWin{i} = Spikes(:,CaTime-30:CaTime+30);
    corr{i} = correlation_dice(syncWin{i});
    figure(i),htmp(corr{i});caxis([0.0 .5]);
end
%total pairwise Synchronization
velCorr = correlation_dice(horzcat(syncWin{:}));
figure,htmp(velCorr);caxis([0.0 .5]);
% Connectivity analysis during locomotion
%pairwise synchornization during rest
threshRest = find(Velocity(:,2)<1);
threshRest = threshRest*4; %multiply by bin width of the encoder time
% Now match threshold to the calcium time by window
for i = 1:length(threshRest)
    CaTime = find(round(time)==threshRest(i));
    syncWinRest{i} = Spikes(:,CaTime-30:CaTime+30);
end
%total pairwise Synchronization
velCorrRest = correlation_dice(horzcat(syncWinRest{:}));
figure,htmp(velCorrRest);caxis([0.0 .5]);
% Connectivity analyis at rest

%% Ensemble Analysis
figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROI,ROIcentroid);
% Ensemble = ensembleNetworks(Ensemble);
%% Plot Ensemble
% ensembleVid(Ensemble,AverageImage,ROIcentroid,files);
% Displays ensemble overlay
[~,I] = sort(cellfun(@length,Ensemble.NodeList),'descend'); %sort max node size
rankEdges = Ensemble.NumEdges(:,I);
rankEnsembles = Ensemble.NodeList(:,I); 
[grad,~]=colorGradient([1 0 0],[0 0 0],4)
figure,
for i = 1:3
    axis off
    color = jet(3);
    EnsembleMap(AverageImage,ROIcentroid,rankEnsembles{i},4,grad(i,:))
    set(gcf,'Position',[100 100 500 500])
    drawnow
    hold on
end
% Combine Maps
figure,imagesc(interp2(Ensemble.sim_index,2)),colormap(jet),caxis([0.13 .4])
K = (1/25)*ones(5);
figure,imagesc(interp2(conv2(Ensemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .3])

%%
sizeE = cellfun(@size,rankEnsembles,'UniformOutput',false);
sizeE = cell2mat(sizeE);
sizeE(sizeE==1) = [];
sizeE = sizeE/max(sizeE);
figure,plot(sizeE)

sizeEdge = cell2mat(rankEdges);
sizeEdge = sizeEdge/max(sizeEdge);
figure,plot(sizeEdge)
%% Ensemble Stats
EnsembleStats
%% Information Entropy
informationEntropy = shannonEntropy(rankEnsembles);
%% SVD/PCA of Ensembles
[tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(Ensemble.sim_index,10);
%% Trial by Trial analysis ##Only use with batch processed files##
addpath(genpath('Figures'));
[batchSpikes,batch_corr] = TrialByTrial(batchData([1,2,4])); % Function call
bin = 20;
[vectorized,sim_index] = cosine_similarity(batchSpikes,bin);
[z,mu,sigma] = zscore(sim_index);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index,100);
caxis([0 0.7]);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
 figure('Name','Dice Correlation')
 for i = 1:size(batch_corr,3)
     subplot(2,3,i),h = htmp(batch_corr(:,:,i),20);caxis([0 0.4]);
 end
%% Plot all the Figures
addpath('Figures');
figure('Name','DeltaF/F'); stack_plot(DeltaFoverF(5:20,1:8000),1.5,8); 
figure('Name','Convolved Spikes'); plot(dDeltaFoverF');
figure('Name','Threshold Detection');DeltaFoverFplotter(dDeltaFoverF,std_threshold,static_threshold)
figure('Name','Spike Plot'); Show_Spikes(Spikes);
% figure('Name','Temporal Shuffled Spike Plot'); shuffledTspikePlot = Show_Spikes(Spikes_shuffled);
% figure('Name','Event Shuffled Spike Plot'); shuffledEspikePlot = Show_Spikes(Event_shuffled);
% figure('Name','Total Shuffled Spike Plot'); shuffledAspikePlot = Show_Spikes(Total_shuffled);
figure('Name','Fluorescence Map'); spike_map(DeltaFoverF);caxis([0 1]),set(gcf,'Position',[100 100 400 400])
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Coactivity Index'); B = bar(coactive_cells,4);ax = gca;ax.TickDir = 'out';ax.Box = 'off';
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.8]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
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
A    = M2;
imwrite(A(:, :, 1), 'test.tiff');
for k = 2:size(A, 3)
  imwrite(A(:, :, k), 'test.tiff', 'WriteMode', 'append');
end
image_movie = mat2gray(M2);
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

