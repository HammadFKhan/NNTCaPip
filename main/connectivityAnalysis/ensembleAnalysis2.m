function Ensemble = ensembleAnalysis2(Spikes,ROIcentroid)
surrogate = 1000;
% Calculate threshold for coactive activations
disp(['Shuffling data ' num2str(surrogate) ' times to find optimal ensemble cutoff']);
shufSpikes = tempShuffle(Spikes,surrogate);
[coactive_cells,~] = coactive_index(Spikes,length(Spikes)/5);
[shufcoactive_cells,~] = coactive_index(shufSpikes,length(Spikes)/5);
bin = ceil(max(coactive_cells)*100);
figure,hold on,bar(coactive_cells),bar(shufcoactive_cells);
% figure,hold on,histogram(coactive_cells,100),histogram(shufcoactive_cells,100);
ensembleWin = 5;
ensembleCan = ensembleWin*find(coactive_cells>(3*std(shufcoactive_cells)+mean(shufcoactive_cells))); % 99% distribution threshold
disp(['Ensemble Candidates: ' num2str(length(ensembleCan))])
% Grab window around ensembles
ensemble = [];
ensembleFrame = []; % reference for what frames were merged for ensemble analysis
for i = 1:length(ensembleCan)
    % Condition if window is outside array bound
    if ensembleCan(i)<=ensembleWin
        ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i):ensembleCan(i)+ensembleWin));
        ensembleFrame = horzcat(ensembleFrame,ensembleCan(i):ensembleCan(i)+ensembleWin);
    elseif size(Spikes,2)-ensembleCan(i)<ensembleWin
        ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i)-ensembleWin:end));
        ensembleFrame = horzcat(ensembleFrame,ensembleCan(i)-ensembleWin:Spikes(end));
        
    else
        ensemble = horzcat(ensemble,Spikes(:,ensembleCan(i)-ensembleWin:ensembleCan(i)+ensembleWin));
        ensembleFrame = horzcat(ensembleFrame,ensembleCan(i)-ensembleWin:ensembleCan(i)+ensembleWin);
    end
end
checkPadding = size(ensemble,1)-size(ensemble,2);
factorCorrection = ensembleWin*floor(size(ensemble,2)/ensembleWin); % Correct for frame size aquisition
fprintf('Vectorizing ensembles candidates...')
if checkPadding>0
    [vectorized,sim_index] = cosine_similarity(horzcat(ensemble,zeros(size(ensemble,1),checkPadding)),1);
    sim_index = sim_index(1:size(ensembleCan,2),1:size(ensembleCan,2));
else
    [vectorized,sim_index] = cosine_similarity(ensemble(:,1:factorCorrection),5);
end
fprintf('done\n')
avgSim = mean(sim_index,'all');
X = sim_index - avgSim*ones(1,size(sim_index,2));
[U,S,V] = svd(X,'econ');

eSVD = cumsum(diag(S))/sum(diag(S)); %energy of svd (find 80% cutoff)
r = 1:find(eSVD>=0.8,1);
Xapprox = U(:,r)*S(r,r)*V(:,r)';
disp(['Sim Index explained by: ' num2str((length(r)/length(eSVD))*100) ' of neural activity'])
figure,
subplot(131),imagesc(sim_index),colormap(jet),caxis([0 1])
subplot(132),imagesc(X),colormap(jet),caxis([0 1])
subplot(133),imagesc(Xapprox),colormap(jet),caxis([0.0 1])
% Weighted Ensemble
ensembleSimilarity = zeros(1,size(S,1));
for i = 1:size(S,1)
    blah = U(:,i)*S(i,i)*V(:,i)';
    ensembleSimilarity(i) = mean(blah(blah>0))+avgSim;
end
% Singular Values
figure,subplot(1,3,1)
semilogx(ensembleSimilarity,'k','Linewidth',2),grid on;hold on;
xlabel('r'); ylabel('Ensemble Similarity')
set(gca,'FontSize',14);
subplot(1,3,2)
semilogx(diag(S),'k','Linewidth',2),grid on;hold on;
xlabel('r'); ylabel('Singular Value, \Sigma_r')
set(gca,'FontSize',14);
subplot(1,3,3)
plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2),grid on;
xlabel('r'); ylabel('Cumulative Energy');
set(gca,'FontSize',14);

% Generate ROC curve using our new Xapprox
thresh = 1;
count = 1;
while thresh>.25 % Cycles threshold value
    xdata(count) = thresh;
    [r,~] = find(tril(Xapprox>thresh,-1));
    rocEnsembles(count) = length(r)+1;
    count = count+1;
    thresh = thresh-0.025;
end
figure,plot(xdata,rocEnsembles,'LineWidth',2)
% Refine ensemble nodes based on similarity and ROC
thresh = .7;
idx = find(round(xdata,2)==thresh);
while rocEnsembles(idx)>100
    idx = idx-1;
    if idx == 0 % exception handling if threshold is never reached
        idx = 1;
        break;
    end
end
thresh = xdata(idx);
disp(['Restricting ensembles to ' num2str(thresh) ' similarity']);
[r,~] = find(tril(Xapprox>thresh,-1));
count = 1;
while isempty(r) || length(r)<2 % Checks to see if we found any ensembles
    thresh = thresh-0.05;
    [r,~] = find(tril(Xapprox>thresh,-1));
    if count>1000 %timer condition if loop becomes infinite
        break;
    end
    count = count+1;
end
if thresh == 0.25 % returns function if thereshold is too low
    disp(['Weak ensembles detected at ' num2str(thresh*100) '% threshold']);
    warning('High false positives detected, stopping analysis...')
    Ensemble = [];
    return;
end
disp(['Ensembles detected at ' num2str(thresh*100) '% threshold']);


r = unique(r);
fEnsemble = find(diff(r)~=1)+1; % Find ensemble index location
if isempty(fEnsemble) % means there is only 1 ensemble with contineous index
    fEnsemble = length(r)+1; % 1 ensemble that ends at the end of index r
    ensembleIdentified = 1;
else
    ensembleIdentified = 1 + length(fEnsemble);
end
fprintf('Combining ensemble trains...\n')
count = 1;
ensembleStability = {};
for i = 1:ensembleIdentified
    % Once ensemble periods are detected find nodes
    if i == 1 
        ensembleId = ensemble(:,r(1:fEnsemble(i)-1)); % Idx ensemble position for the first seperation
    elseif i == ensembleIdentified && ensembleIdentified>1
        ensembleId = ensemble(:,r(fEnsemble(i-1):end)); % Idx ensemble position for last seperation
    else
        ensembleId = ensemble(:,r(fEnsemble(i-1):fEnsemble(i)-1)); % Idx ensemble position for all other seperation
    end
    fprintf(['Connectivity analysis for activation ' num2str(i) ' of ' num2str(ensembleIdentified) '.\n']);
    corr = correlation_dice(ensembleId);
    thresh = 0.5;
    Connected_ROI{i} = Connectivity_dice(corr,thresh);
    [NumActiveNodes,NodeList{i},NumNodes{i},NumEdges{i},SpatialCentroid{i},SpatialCentroidVariance{i},...
        ActivityCentroid{i},ActivityCentroidVariance{i}, ActivityCoords{i}]...
        = Network_Analysis(ROIcentroid,Connected_ROI{i});
    if ~isempty(NodeList{i})
        try
            ensembleStability{count} = coactive_index(ensemble(NodeList{i},:),length(ensemble));
            count = count+1;
        catch ME
            continue;
        end
    end
end
% % Norm ensemble stability
% for i = 1:size(ensembleStability,2)
%     ensembleStability(:,i) = ensembleStability(:,i)/(max(ensembleStability(:,i))-min(ensembleStability(:,i)));
% end
% ensembleStability = [];
fprintf('done\n')
Ensemble.ensemble = ensemble;
Ensemble.ensembleCan = ensembleCan;
Ensemble.ensembleFrame = ensembleFrame;
Ensemble.ensembleStability = ensembleStability;
Ensemble.vectorized = vectorized;
Ensemble.sim_index = sim_index;
Ensemble.sim_indexApprox = Xapprox;
% Ensemble.NetworkAnalysis = NetworkAnalysis;
Ensemble.NumActiveNodes = NumActiveNodes;
Ensemble.NodeList = NodeList;
Ensemble.NumNodes = NumNodes;
Ensemble.NumEdges = NumEdges;
Ensemble.SpatialCentroid = SpatialCentroid;
Ensemble.SpatialCentroidVariance = SpatialCentroidVariance;
Ensemble.ActivityCentroid = ActivityCentroid;
Ensemble.ActivityCentroidVariance = ActivityCentroidVariance;
Ensemble.ActivityCoords =  ActivityCoords;
Ensemble.ensembleIndentified = ensembleIdentified;
Ensemble.Connected_ROI =  Connected_ROI;
Ensemble.loc = r;
Ensemble.coactive_cells = coactive_cells;
Ensemble.shufcoactive_cells = shufcoactive_cells;
disp(['Ensembles Identified: ' num2str(ensembleIdentified)])
end



