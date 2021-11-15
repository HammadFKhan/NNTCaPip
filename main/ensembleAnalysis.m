function Ensemble = ensembleAnalysis(Spikes,ROI,ROIcentroid)
surrogate = 10000;
% Calculate threshold for coactive activations
disp(['Shuffling data ' num2str(surrogate) ' times to find optimal ensemble cutoff']);
shufSpikes = tempShuffle(Spikes,200);
[coactive_cells,~] = coactive_index(Spikes,length(Spikes));
[shufcoactive_cells,~] = coactive_index(shufSpikes,length(Spikes));
bin = ceil(max(coactive_cells)*100);
% hold on,bar(coactive_cells),bar(shufcoactive_cells);
ensembleCan = find(coactive_cells>(2.5*std(shufcoactive_cells)+mean(shufcoactive_cells))); % 99% distribution threshold
for i = 1:length(ensembleCan)
    ensemble(:,i) = Spikes(:,ensembleCan(i));
end
checkPadding = size(ensemble,1)-size(ensemble,2);
if checkPadding>0
    [vectorized,sim_index] = cosine_similarity(horzcat(ensemble,zeros(size(ensemble,1),checkPadding)),1);
    sim_index = sim_index(1:size(ensembleCan,2),1:size(ensembleCan,2));
else
    [vectorized,sim_index] = cosine_similarity(ensemble,5);
end

% Refine ensemble nodes based on similarity
thres = 0.6;
[r,~] = find(tril(sim_index>thres,-1));
if isempty(r)
    disp(['No ensembles detected at ' num2str(thres*100) '% threshold']);
    return;
end

r = unique(r);
fEnsemble = find(diff(r)~=1)+1;
ensembleIndentified = 1 + length(fEnsemble);

for i = 1:ensembleIndentified
    % Once ensemble periods are detected find nodes
    if i == 1
        ensembleId = ensemble(:,r(1:fEnsemble(i)-1)); % Idx ensemble position for the first seperation
    elseif i == ensembleIndentified
        ensembleId = ensemble(:,r(fEnsemble(i-1):end)); % Idx ensemble position for last seperation
    else
        ensembleId = ensemble(:,r(fEnsemble(i-1):fEnsemble(i)-1)); % Idx ensemble position for all other seperation
    end
    
    corr = correlation_dice(ensembleId);
    Connected_ROI{i} = Connectivity_dice(corr, ROI);
    [NumActiveNodes,NodeList{i},NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
        ActivityCentroid,ActivityCentroidVariance]...
        = Network_Analysis(ROIcentroid,Connected_ROI{i});
end
Ensemble.ensemble = ensemble;
Ensemble.sim_index = sim_index;
Ensemble.NumActiveNodes = NumActiveNodes;
Ensemble.NodeList = NodeList;
Ensemble.NumNodes = NumNodes;
Ensemble.NumEdges = NumEdges;
Ensemble.SpatialCentroid = SpatialCentroid;
Ensemble.SpatialCentroidVariance = SpatialCentroidVariance;
Ensemble.ActivityCentroid = ActivityCentroid;
Ensemble.ActivityCentroidVariance = ActivityCentroidVariance;
Ensemble.ensembleIndentified = ensembleIndentified;
Ensemble.loc = r;
disp(['Ensembles Identified: ' num2str(ensembleIndentified)])
end



