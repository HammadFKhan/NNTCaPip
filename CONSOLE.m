%% NNT Calcium Pipeline
% Github Version 4.1
%% Remove ROIs
if exist('badComponents','var') && ~exist('badComFlag','var')
    [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
        removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
    badComFlag = 1;
end
%% Fix centroids
ROIcentroid = [];
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    ROIcentroid(i,:) = floor(mean(blah,1));
end
%% Analysis
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'));
std_threshold = 2;
static_threshold = 0;
Spikes = Spike_Detector_Single(diff(DeltaFoverF),std_threshold,static_threshold);
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,size(Spikes,2));
cell_count = length(ROI);
time = time_adjust(size(DeltaFoverF,2),30.048);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,250);%std, window (frames)
end

% Perform shuffling and pairwise if data is small enough
if size(DeltaFoverF,2)<2000    
    Spikes_shuffled = tempShuffle(Spikes,1000);
    Event_shuffled = spatialShuffle(Spikes,1000);
    surrogate = 10;
    Total_shuffled = allShuffle(Spikes,1000);
    [shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
    shuff_corr = correlation_dice(Event_shuffled);
    [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
    shufsim_index = shufsim_index-mean(mean(shufsim_index,2));
    factorCorrection = 10*floor(size(Spikes,2)/10); % Correct for frame size aquisition
    [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),10);
    corr = correlation_dice(Spikes);
    Connected_ROI = Connectivity_dice(corr,0.1);
    [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
        ActivityCentroid,ActivityCentroidVariance]...
        = Network_Analysis(ROIcentroid,Connected_ROI);
end
% Pairwise Velocity Analysis
% velocityPairwise(VR_data,Spikes)

%% Ensemble Analysis
% figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);

% Ensemble stats
Ensemble = ensembleMetric(Ensemble,AverageImage,ROIcentroid);
Ensemble = ensembleStat(Ensemble);
%%
W12_10Entropy.informationEntropy = informationEntropy;
W12_10Entropy.rankedEnsembles = rankEnsembles;
W12_10Entropy.rankedEdges = rankEdges;
W12_10Entropy.Ensemble = Ensemble;
%% LFP pipette analysis
addpath(genpath('C:\Users\khan332\Documents\GitHub\NNTEphysPip\main'))
LFP = Ca_LFP(time,1); %caTime; loadFlag0/1; LFP.out
%% Beta Analysis
[peakAlign,norm,f,stats] = IntrabetaAnalysis(LFP.beta);
figure
for i = 1:144
    subplot(12,12,i),plot(LFP.beta.betaTrace{i}),axis off
end

figure,plot(0:1/LFP.Fs:(length(LFP.Vmfilt)-1)/LFP.Fs,LFP.Vmfilt);xlim([0 length(LFP.Vmfilt)/LFP.Fs])

figure,plot(0:1/LFP.Fs:(length(LFP.betaLFP)-1)/LFP.Fs,LFP.betaLFP);xlim([0 length(LFP.betaLFP)/LFP.Fs])

%% Behavior
Vel = EncoderVelocity2(abs(encoder(:,1)),abs(encoder(:,2))); % position;time
%% Generate Rest/Run Ca Spikes
thresh = 2;
if ~exist('CaFR','var'), CaFR = 30.048;end % sets to default framerate
if ~exist('time','var'), timeT = 1/CaFR:1/CaFR:(size(DeltaFoverF,2)/CaFR);end
[runSpikes,runSpikesFrame] = spikeState(Vel,Spikes,time,CaFR,thresh,1); % state 1/0 for run/rest
[restSpikes,restSpikesFrame] = spikeState(Vel,Spikes,time,CaFR,thresh,0); % state 1/0 for run/rest

if iscell(runSpikes)
    runSpikes = horzcat(runSpikes{:});
end
if iscell(restSpikes)
   restSpikes = horzcat(restSpikes{:});
end
runCa = DeltaFoverF(:,runSpikesFrame);
restCa = DeltaFoverF(:,restSpikesFrame);

%% Behavior-based Ensemble
runfactorCorrection = 40*floor(size(runSpikes,2)/40); % Correct for frame size aquisition
restfactorCorrection = 20*floor(size(restSpikes,2)/20); % Correct for frame size aquisition

[vectorizedRun,sim_indexRun] = cosine_similarity(runSpikes(:,1:runfactorCorrection),40);
[vectorizedRest,sim_indexRest] = cosine_similarity(restSpikes(:,1:restfactorCorrection),20);
%%
runEnsemble = ensembleAnalysis(runSpikes(:,1:runfactorCorrection),ROI,ROIcentroid);
restEnsemble = ensembleAnalysis(restSpikes(:,1:restfactorCorrection),ROI,ROIcentroid);

ensembleMetric(runEnsemble,AverageImage,ROIcentroid)
ensembleMetric(restEnsemble,AverageImage,ROIcentroid)


%% SVD/PCA of Ensembles
comVect = [runCa restCa];
[X] = featureProject(comVect);
%% K-means clustering of neural projections
[idx,C] = kmeans(X,2);
scatter3(X(idx==1,1),X(idx==1,2),X(idx==1,3),10,[0 0 0],'filled'); hold on; %[43 57 144]/255
scatter3(X(idx==2,1),X(idx==2,2),X(idx==2,3),10,[1 0 0],'filled'); %[0 148 68]/255
%scatter3(X(idx==3,1),X(idx==3,2),X(idx==3,3),10,[0 0 0],'filled'); %[0 148 68]/255

%% Beta events and fluorescence


%% Beta events within ensembles
runLFP = betaCaEnsemble(runSpikes,runSpikesFrame,runEnsemble,LFP,CaFR); 

restLFP = betaCaEnsemble(restSpikes,restSpikesFrame,restEnsemble,LFP,CaFR); 

[peakAlign,norm,f,stats] = IntrabetaAnalysis(runLFP.beta);
[peakAlign,norm,f,stats] = IntrabetaAnalysis(restLFP.beta);

%% Plot all the Figures
addpath('Figures');
figure('Name','DeltaF/F'); stack_plot(DeltaFoverF,1.5,15); 
figure('Name','Convolved Spikes'); plot(dDeltaFoverF');
figure('Name','Threshold Detection');DeltaFoverFplotter(dDeltaFoverF,std_threshold,static_threshold)
figure('Name','Spike Plot'); Show_Spikes(Spikes);
% figure('Name','Temporal Shuffled Spike Plot'); shuffledTspikePlot = Show_Spikes(Spikes_shuffled);
% figure('Name','Event Shuffled Spike Plot'); shuffledEspikePlot = Show_Spikes(Event_shuffled);
% figure('Name','Total Shuffled Spike Plot'); shuffledAspikePlot = Show_Spikes(Total_shuffled);
figure('Name','Fluorescence Map'); spike_map(DeltaFoverF);caxis([0 1])
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);
figure('Name','Coactivity Index'); B = bar(coactive_cells,4);ax = gca;ax.TickDir = 'out';ax.Box = 'off';
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Dice-Similarity Index');h = htmp(shuff_corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.0 .8]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Cosine-Similarity Index'); h = htmp(shufsim_index);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Network Map'); NodeSize = 5;EdgeSize = 1;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,NodeSize,EdgeSize)


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
[tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(vectorized,139);
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
 
plot_raster(1:120,Spikes(5,1:120))
% Have you tried using Multidimensional Scaling (MDS) to emebed the
% centroids in a 2 dimensional space for visualization?

% This should visualize how the centroids related to each other. You could�
% also then compute the Delauney Triangulation of the projected graph, to
% identify neighbors.

