%% NNT Calcium Pipeline
% Github Version 4.1
%% Remove ROIs
if exist('badComponents','var')
[DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
    removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
end
%% Analysis
addpath(genpath('main'));
std_threshold = 8;
static_threshold = .01;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,500);
cell_count = length(ROI);
time = time_adjust(size(DeltaFoverF,2),30.048);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,100);%std, window (frames)
end

% Perform shuffling and pairwise if data is small enough
if size(DeltaFoverF,2)<2000    
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
% Pairwise Velocity Analysis
% velocityPairwise(VR_data,Spikes)
% Ensemble Analysis
% figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROI,ROIcentroid);
% Ensemble = ensembleNetworks(Ensemble);
% Plot Ensemble
% ensembleVid(Ensemble,AverageImage,ROIcentroid,files);
% Displays ensemble overlay
[~,I] = sort(cellfun(@length,Ensemble.NodeList),'descend'); %sort max node size
rankEdges = Ensemble.NumEdges(:,I);
rankEnsembles = Ensemble.NodeList(:,I); 
[grad,~]=colorGradient([1 0 0] ,[0 0 0],5)
Ensemble.rankEnsembles = rankEnsembles;
figure,
for i = 1:3
    axis off
    color = jet(3);
    EnsembleMap(AverageImage,ROIcentroid,rankEnsembles{i},5,grad(i,:))
    set(gcf,'Position',[100 100 500 500])
    drawnow
    hold on
end
% Combine Maps
figure,imagesc(interp2(Ensemble.sim_index,2)),colormap(jet),caxis([0.13 .4])
K = (1/25)*ones(5);
figure,imagesc(interp2(conv2(Ensemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .3])

%
sizeE = cellfun(@size,rankEnsembles,'UniformOutput',false);
sizeE = cell2mat(sizeE);
sizeE(sizeE==1) = [];
sizeE = sizeE/max(sizeE);
figure,plot(sizeE),title('Ranked Ensembles')

sizeEdge = cell2mat(rankEdges);
sizeEdge = sizeEdge/max(sizeEdge);
figure,plot(sizeEdge),title('Ranked Connections')
% Ensemble Stats
% EnsembleStats
% Information Entropy
informationEntropy = shannonEntropy(rankEnsembles);
% Plot Centroid Boundary
figure,hold on
%Rank by number of Activity points
[~,I] = sort(cellfun(@length,Ensemble.ActivityCoords),'descend'); %sort max node size
rankedActivityCoords = Ensemble.ActivityCoords(:,I);
Ensemble.rankedActivityCoords = rankedActivityCoords;
checkSize = cell2mat(cellfun(@size,rankedActivityCoords,'UniformOutput',false)');
for i = 1:size(checkSize(checkSize(:,1)>2),1)
    x = rankedActivityCoords{i}(:,1);y = rankedActivityCoords{i}(:,2);
    k = boundary(x,y);
    x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
    y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
    x2 = smoothdata(x1,'gaussian',50);
    y2 = smoothdata(y1,'gaussian',50);
    plot(x2,y2,'Color',[0.5 0.5 0.5])
    scatter(x,y,4,'k','filled')
end

%%
W12_10Entropy.informationEntropy = informationEntropy;
W12_10Entropy.rankedEnsembles = rankEnsembles;
W12_10Entropy.rankedEdges = rankEdges;
W12_10Entropy.Ensemble = Ensemble;
%% LFP pipette analysis
addpath('C:\Users\khan332\Documents\GitHub\NNTEphysPip')
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
Vel = EncoderVelocity2(encoder(:,1),encoder(:,2)); % position;time
%%Generate Rest/Run Ca Spikes
[runSpikes,runSpikesFrame] = spikeState(Vel,Spikes,time,CaFR,1); % state 1/0 for run/rest
[restSpikes,restSpikesFrame] = spikeState(Vel,Spikes,time,CaFR,0); % state 1/0 for run/rest

if iscell(runSpikes)
    runSpikes = horzcat(runSpikes{:});
end
if iscell(restSpikes)
   restSpikes = horzcat(restSpikes{:});
end


%% Behavior-based Ensemble
runfactorCorrection = 5*floor(size(runSpikes,2)/5); % Correct for frame size aquisition
restfactorCorrection = 5*floor(size(restSpikes,2)/5); % Correct for frame size aquisition

[vectorizedRun,sim_indexRun] = cosine_similarity(runSpikes(:,1:runfactorCorrection),5);
[vectorizedRest,sim_indexRest] = cosine_similarity(restSpikes(:,1:restfactorCorrection),25);

runEnsemble = ensembleAnalysis(runSpikes(:,1:runfactorCorrection),ROI,ROIcentroid);
restEnsemble = ensembleAnalysis(restSpikes(:,1:restfactorCorrection),ROI,ROIcentroid);

ensembleMetric(runEnsemble,AverageImage,ROIcentroid)
ensembleMetric(restEnsemble,AverageImage,ROIcentroid)


%% SVD/PCA of Ensembles
comVect = [vectorizedRun vectorizedRest];
[tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(comVect,length(vectorizedRun));
%%
runLFP = betaCaEnsemble(runSpikes,runSpikesFrame,runEnsemble,LFP,CaFR); 


restLFP = betaCaEnsemble(restSpikes,restSpikesFrame,restEnsemble,LFP,CaFR); 
%% Beta events within ensembles
frameT = runSpikesFrame(runEnsemble.ensembleFrame);
count = 1;
for i = 1:286
    temp = find(frameT(i)==betaEventFrame(:,3)); % checks to see if a beta event lies on the frame\
    if temp
        disp(['Beta event matched to idx: ' num2str(i)])
        hold on,xline(i,'r');
        LFP.beta.ensembleBetaMatch(count,1) = frameT(i); %frame
        LFP.beta.ensembleBetaMatch(count,2) = temp; %beta index
        count = count+1;
    end
end

runLFP.beta = LFP.beta; % create a structure array looking at only behavior based beta
runLFP.beta.betaBurst.detectedBeta = LFP.beta.betaBurst.detectedBeta(LFP.beta.ensembleBetaMatch(:,2),:);

% plot beta traces
figure,
for i = 1:30
    subplot(5,6,i),plot(LFP.beta.betaTrace{restLFP.beta.ensembleBetaMatch(i,2)}), axis off
end

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
figure('Name','Fluorescence Map'); spike_map(fluorescenceTrials(:,:,3));caxis([0 1]),set(gcf,'Position',[100 100 400 400])
figure('Name','Population Intensity');height = 10;rateImage = firing_rate(Spikes,height,time);caxis([0 0.5]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Coactivity Index'); B = bar(coactive_cells,4);ax = gca;ax.TickDir = 'out';ax.Box = 'off';
figure('Name','Dice-Similarity Index');h = htmp(corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Dice-Similarity Index');h = htmp(shuff_corr,10);caxis([0 0.4]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Cosine-Similarity Index'); h = htmp(sim_index);caxis([0.35 .9]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Shuffled Cosine-Similarity Index'); h = htmp(shufsim_index);caxis([0 1]);set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
figure('Name','Network Map'); NodeSize = 0;EdgeSize = 1;Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid,ROI,NodeSize,EdgeSize)


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

