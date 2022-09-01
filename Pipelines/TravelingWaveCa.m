%Traveling Waves Planer Grids + Imaging

ROIcentroid = {json_file.centroid}';ROIcentroid = cell2mat(ROIcentroid);
for i = 1:size(json_file,1)
    for ii = 1:length(json_file(i).coordinates)
        ROI{i,1}{ii,1} = json_file(i).coordinates(ii,1:2)';
    end
end
DeltaFoverF = full(C_df);
idx = size(C_df,1);
parfor i = 1:idx
    disp('Deconvolving using contrained metric...')
    [c(i,:), s(i,:), options] = deconvolveCa(full(C_df(i,:)));
end
dDeltaFoverF = s;
%% Analysis
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
time = time_adjust(size(DeltaFoverF,2),30);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,120);
end

% Perform shuffling and pairwise if data is small enough
if size(DeltaFoverF,2)<2000    
%     Spikes_shuffled = tempShuffle(Spikes,1000);
%     Event_shuffled = spatialShuffle(Spikes,1000);
%     surrogate = 10;
%     Total_shuffled = allShuffle(Spikes,1000);
%     [shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
%     shuff_corr = correlation_dice(Event_shuffled);
%     [shufvectorized,shufsim_index] = cosine_similarity(Spikes_shuffled,10);
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
[grad,~]=colorGradient([1 0 0],[0 0 0],3);
Ensemble.rankEnsembles = rankEnsembles;
figure,
for i = 1
    axis off
    color = jet(3);
    EnsembleMap(AverageImage,ROIcentroid,lateSpikeEnsemble.rankEnsembles{i},4,grad(i,:))
    set(gcf,'Position',[100 100 500 500])
    drawnow
    hold on
end
% Combine Maps
figure,imagesc(interp2(Ensemble.sim_index,2)),colormap(jet),caxis([0.13 .4])
K = (1/25)*ones(5);
figure,imagesc(interp2(conv2(Ensemble.sim_index,K,'same'),2)),colormap(jet),caxis([.08 .3])

%
sizeE = cellfun(@size,Ensemble.rankEnsembles,'UniformOutput',false);
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
informationEntropy = shannonEntropy(Ensemble.rankEnsembles);
% Plot Centroid Boundary
figure,hold on
%Rank by number of Activity points
[~,I] = sort(cellfun(@length,Ensemble.ActivityCoords),'descend'); %sort max node size
rankedActivityCoords = Ensemble.ActivityCoords(:,I);
Ensemble.rankedActivityCoords = rankedActivityCoords;
% for i = 1:15%size(Ensemble.ActivityCoords,2)
%     
%     x = rankedActivityCoords{i}(:,1);y = rankedActivityCoords{i}(:,2);
%     k = boundary(x,y);
%     x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
%     y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
%     x2 = smoothdata(x1,'gaussian',50);
%     y2 = smoothdata(y1,'gaussian',50);
%     plot(x2,y2,'Color',[0.5 0.5 0.5])
%     scatter(x,y,4,'k','filled')
%     
% end
%%
ensembleVid(nolateSpikeEnsemble,AverageImage,ROIcentroid,'nolateSpikeEnsemble')
%% Now look at trial specific changes
% Standard trial length is 209 frames
trialLength = 360;
for i = 1:floor(size(Spikes,2)/trialLength)
    spikeTrials{i} = Spikes(:,((i-1)*trialLength+1):i*trialLength);
end
    
% figure,
% for i = 1:length(spikeTrials)
%     subplot(6,6,i),Show_Spikes(spikeTrials{i})
% end


figure,
for i = 1:floor(size(Spikes,2)/trialLength)
    Show_Spikes(spikeTrials{i}),hold on
end
%% generate ROC for connectivity
thresh = 0;
connectedNum = [];
step = 0.025;
for i = 1:(1/step)+1
    Connected_ROI = Connectivity_dice(corr,ROI,thresh);
    connectedNum(i) = size(Connected_ROI,1);
    thresh = thresh+step;
end
figure,plot(0:step:1,connectedNum);
%% Trial specific connectivity
simM = [];
Connected_ROI = [];
Ls = [];
figure,
count = 1;
for i = trialData.responsiveTrials.lateSpikeTrials  %t1
    corr = correlation_dice(spikeTrials{i});
    Connected_ROI{i} = Connectivity_dice(corr,ROI,0.4);
    [NumActiveNodes{i},NodeList{i},NumNodes{i},NumEdges{i},SpatialCentroid{i},SpatialCentroidVariance{i},ActivityCentroid{i},ActivityCentroidVariance{i}, ActivityCoords{i}]...
        = Network_Analysis(ROIcentroid,Connected_ROI{i});
    [vectorized,sim_index] = cosine_similarity(spikeTrials{i},10);
    subplot(5,4,count),imagesc(sim_index),colormap(jet)
    Ls(count) = mean(sim_index(sim_index>0.1),'all'); %mean center
    count = count+1;
end
% t1 = vertcat(simValue{:});
 figure,boxplot(Ls);ylim([0.0 0.5])
%%
Connected_ROI = {};
simM = [];
simValue = [];
nLs = [];
count = 1;
figure,
for i = trialData.responsiveTrials.noLateSpikeTrials%t2
    corr = correlation_dice(spikeTrials{i});
    Connected_ROI{i} = Connectivity_dice(corr,ROI,0.3);
    [NumActiveNodes{i},NodeList{i},NumNodes{i},NumEdges{i},SpatialCentroid{i},SpatialCentroidVariance{i},ActivityCentroid{i},ActivityCentroidVariance{i}, ActivityCoords{i}]...
        = Network_Analysis(ROIcentroid,Connected_ROI{i});
    [vectorized,sim_index] = cosine_similarity(spikeTrials{i},10);
    subplot(5,4,count),imagesc(sim_index),colormap(jet)
    nLs(count) = mean(sim_index(sim_index>0),'all');
    count = count+1;
end
% t2 = vertcat(simValue{:});
figure,boxplot(nLs);ylim([0.0 0.5])
%%
[vectorized,sim_index] = cosine_similarity(horzcat(spikeTrials{trialData.responsiveTrials.lateSpikeTrials(1:5)} ),10);
%%
count = 1;
for i = 1:49
    if ~isempty(lateSpikeEnsemble.ActivityCentroid{i}) 
        euclideanCentroid(count) = sqrt(lateSpikeEnsemble.ActivityCentroid{i}(1)^2+lateSpikeEnsemble.ActivityCentroid{i}(2)^2);
        count = count+1;
    end
end
%%
count = 1;
for i = 1:49
    if ~isempty(nolateSpikeEnsemble.ActivityCentroid{i}) 
        euclideanCentroid(count) = sqrt(nolateSpikeEnsemble.ActivityCentroid{i}(1)^2+nolateSpikeEnsemble.ActivityCentroid{i}(2)^2);
        count = count+1;
    end
end
%%
NodeSize = 2;EdgeSize = 1;
figure('Name','Network Map'),hold on
for i = 1:length(spikeTrials)
    subplot(7,7,i),htmp(corr{26},10),title(['Trial ' num2str(i)])
end
%% Quantify Node Reactivation
Connected_ROI = []
for i = trialData.responsiveTrials.noLateSpikeTrials%t2
    corr = correlation_dice(spikeTrials{i});
    Connected_ROI{i} = Connectivity_dice(corr,ROI,0.3);
end

nodeWeight = [];
NodeProbability = [];
for i = 1:size(Spikes,1)
    count = 1;
    temp = [];
    for ii = 1:length(Connected_ROI)
        if ~isempty(Connected_ROI{ii})
            [r,~] = find(Connected_ROI{ii}(:,1:2)==i);
            temp(count) = mean(length(r)/size(Connected_ROI{ii},1));
            temp2(count) = mean(Connected_ROI{ii}(r,3));
            if isnan(temp2(count))
                temp2(count) = 0;
            end
            count = count+1;
        end
    end
    NodeWeight(i) = mean(temp2);
    NodeProbability(i) = mean(temp);
end

% Normalize node weights and find threshold for crit nodes
normNodeWeight = normalize(NodeWeight,'range');
normNodeProbability = normalize(NodeProbability,'range');
nodeWThresh = 0.6;
nodePThresh = 0.6;

% find crit nodes across each axis
critWNodes = find(normNodeWeight>nodeWThresh);
critWNodeValue = normNodeWeight(normNodeWeight>nodeWThresh);
critPNodes = find(normNodeProbability>nodePThresh);
critPNodeValue = normNodeProbability(normNodeProbability>nodePThresh);
% check membership for each 
if length(critWNodes)>length(critPNodes) %checks whick one is larger for membership assignment
    critNodes = critPNodes(ismember(critPNodes,critWNodes)); 
else 
    critNodes = critWNodes(ismember(critWNodes,critPNodes));
end

figure,plot(NodeWeight,NodeProbability,'.');
xline(nodeWThresh,'k--');
yline(nodePThresh,'k--');
%% Late vs no Late spike ensembles
[lateSpikeEnsemble, nolateSpikeEnsemble] =...
    travelingWaveEnsemble(spikeTrials,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROI,ROIcentroid,AverageImage);

%% some statistics about these ensembles
lateSpikeEnsemble = ensembleStat(lateSpikeEnsemble);
nolateSpikeEnsemble = ensembleStat(nolateSpikeEnsemble);
%% manifold analysis and entropy
lateSpikeEnsemble = ensembleMetric(lateSpikeEnsemble,AverageImage,ROIcentroid);
nolateSpikeEnsemble = ensembleMetric(nolateSpikeEnsemble,AverageImage,ROIcentroid);
%% Overlapping ensemble metric
% Check number of shared nodes by adding each interative ensemble rank
checkIteration = length(lateSpikeEnsemble.rankEnsembles)- length(nolateSpikeEnsemble.rankEnsembles);
if checkIteration>0 %ie more late spike ensembles
    ensembleIteration = length(nolateSpikeEnsemble.rankEnsembles); %only index across smaller one (and we'll negate the remaining number
else %ie more no late spike ensembles
    ensembleIteration = length(lateSpikeEnsemble.rankEnsembles);
end

sharedEnsemble = [];
count = 1;
for i = 1:ensembleIteration
    lsEnsemble = lateSpikeEnsemble.rankEnsembles{i};
    nlsEnsemble = nolateSpikeEnsemble.rankEnsembles{i};
    checkIdx = length(lsEnsemble)-length(nlsEnsemble); %find which one is larger and check existing indices
    if ~isempty(checkIdx)
        if checkIdx>0
            idx = ismember(lsEnsemble,nlsEnsemble);
        else
            idx = ismember(lsEnsemble,nlsEnsemble);
        end
        if count>1
            sharedEnsemble(count) = sum(idx)+sharedEnsemble(count-1);
        else
            sharedEnsemble(count) = sum(idx);
        end
        count = count+1;
    end
end

%%
lateSpikeConnectedROI = vertcat(lateSpikeEnsemble.Connected_ROI{:});
[r,c] = find(lateSpikeConnectedROI==0);
lateSpikeConnectedROI(r,:) = [];
nolateSpikeConnectedROI = vertcat(nolateSpikeEnsemble.Connected_ROI{:});
[r,c] = find(nolateSpikeConnectedROI==0);
nolateSpikeConnectedROI(r,:) = [];
figure('Name','Late Spike Network Map'); NodeSize = 3;EdgeSize = 1;Cell_Map_Dice(AverageImage,lateSpikeConnectedROI,ROIcentroid,NodeSize,EdgeSize)
figure('Name','No Late Spike Network Map'); NodeSize = 2;EdgeSize = 2;Cell_Map_Dice(AverageImage,nolateSpikeConnectedROI,ROIcentroid,NodeSize,EdgeSize)
%%
W12_10Entropy.informationEntropy = informationEntropy;
W12_10Entropy.rankedEnsembles = rankEnsembles;
W12_10Entropy.rankedEdges = rankEdges;
W12_10Entropy.Ensemble = Ensemble;

%% SVD/PCA of Ensembles
[vectorizedL,sim_indexL] = cosine_similarity(lateSpikeEnsemble.Spikes,10);
[vectorizedNL,sim_indexNL] = cosine_similarity(nolateSpikeEnsemble.Spikes,40);

comVect = [vectorizedL vectorizedNL];
[tProjq1, tProjq2, uProjq1, uProjq2] = featureProject(comVect,length(vectorizedL),0);
legend('Late Spike','No Late Spike')

comSim = [sim_indexL sim_indexL];
svd_analysis(comSim)
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
figure,plot(strongConnections(:,8))
hold on,plot(weakConnections(:,8))
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4 3]);
%%
plot_raster(1:120,Spikes(5,1:120))
% Have you tried using Multidimensional Scaling (MDS) to emebed the
% centroids in a 2 dimensional space for visualization?


