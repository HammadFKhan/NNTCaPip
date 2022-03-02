function Ensemble = ensembleAnalysis(Spikes,ROI,ROIcentroid)
surrogate = 1000;
% Calculate threshold for coactive activations
disp(['Shuffling data ' num2str(surrogate) ' times to find optimal ensemble cutoff']);
shufSpikes = tempShuffle(Spikes,surrogate);
[coactive_cells,~] = coactive_index(Spikes,length(Spikes)/5);
[shufcoactive_cells,~] = coactive_index(shufSpikes,length(Spikes)/5);
bin = ceil(max(coactive_cells)*100);
figure,hold on,bar(coactive_cells),bar(shufcoactive_cells);
% figure,hold on,histogram(coactive_cells,100),histogram(shufcoactive_cells,100);
ensembleCan = 5*find(coactive_cells>(2.5*std(shufcoactive_cells)+mean(shufcoactive_cells))); % 99% distribution threshold
disp(['Ensemble Candidates: ' num2str(length(ensembleCan))])
% Grab window around ensembles
ensembleWin = 5;
ensemble = [];
for i = 1:length(ensembleCan)
    ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i)-ensembleWin:ensembleCan(i)+ensembleWin));
    % Condition if window is outside array bound
    if ensembleCan<ensembleWin
            ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i):ensembleCan(i)+ensembleWin));
    elseif (size(Spikes,2)-ensembleCan)<ensembleWin
            ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i)-ensembleWin:end));
    end
end
checkPadding = size(ensemble,1)-size(ensemble,2);
factorCorrection = ensembleWin*floor(size(ensemble,2)/ensembleWin); % Correct for frame size aquisition
if checkPadding>0
    [vectorized,sim_index] = cosine_similarity(horzcat(ensemble,zeros(size(ensemble,1),checkPadding)),1);
    sim_index = sim_index(1:size(ensembleCan,2),1:size(ensembleCan,2));
else
    [vectorized,sim_index] = cosine_similarity(ensemble(:,1:factorCorrection),ensembleWin);
end

% Generate ROC curve
thresh = 1;
[r,~] = find(tril(sim_index>thresh,-1));
count = 1;
while thresh>.25 % Cycles threshold value
    count = count+1;
    thresh = thresh-0.025;
    [r,~] = find(tril(sim_index>thresh,-1));
    rocEnsembles(count) = length(r)+1;
end
figure,plot(flip(.25:.025:1),rocEnsembles,'LineWidth',2)
% Refine ensemble nodes based on similarity and ROC
thresh = .6;
[r,~] = find(tril(sim_index>thresh,-1));
count = 1;
while isempty(r) || length(r)<2 % Checks to see if we found any ensembles
    thresh = thresh-0.05;
    [r,~] = find(tril(sim_index>thresh,-1));
    if count>1000 %timer condition if loop becomes infinite
        break;
    end
    count = count+1;
end
if thresh == 0.25 % returns function if thereshold is too low
    disp(['Weak ensembles detected at ' num2str(thresh*100) '% threshold']);
    warning('High false positives detected, stopping analysis...')
%     Ensemble = [];
%     return;
end
disp(['Ensembles detected at ' num2str(thresh*100) '% threshold']);

r = unique(r);
fEnsemble = find(diff(r)~=1)+1; % Find ensemble index location
if isempty(fEnsemble) % means there is only 1 ensemble with contineous index
    fEnsemble = length(r)+1; % 1 ensemble that ends at the end of index r
    ensembleIndentified = 1;
else
    ensembleIndentified = 1 + length(fEnsemble);
end

for i = 1:ensembleIndentified
    % Once ensemble periods are detected find nodes
    if i == 1 
        ensembleId = ensemble(:,r(1:fEnsemble(i)-1)); % Idx ensemble position for the first seperation
    elseif i == ensembleIndentified && ensembleIndentified>1
        ensembleId = ensemble(:,r(fEnsemble(i-1):end)); % Idx ensemble position for last seperation
    else
        ensembleId = ensemble(:,r(fEnsemble(i-1):fEnsemble(i)-1)); % Idx ensemble position for all other seperation
    end
    
    corr = correlation_dice(ensembleId);
    thresh = 0.2;
    Connected_ROI{i} = Connectivity_dice(corr, ROI,thresh);
    [NumActiveNodes,NodeList{i},NumNodes{i},NumEdges{i},SpatialCentroid{i},SpatialCentroidVariance{i},...
        ActivityCentroid{i},ActivityCentroidVariance{i}]...
        = Network_Analysis(ROIcentroid,Connected_ROI{i});
end
Ensemble.ensemble = ensemble;
Ensemble.vectorized = vectorized;
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
Ensemble.Connected_ROI =  Connected_ROI;
Ensemble.loc = r;
disp(['Ensembles Identified: ' num2str(ensembleIndentified)])
end



