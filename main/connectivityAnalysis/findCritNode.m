function Ensemble = findCritNode(trials,ROI,spikeTrials,Ensemble)

Connected_ROI = [];
count = 1;
fprintf('Correlating trials...')
for i = trials%t2
    corr = correlation_dice(spikeTrials{i});
    Connected_ROI{count} = Connectivity_dice(corr,ROI,0.3);
    count = count+1;
end
fprintf('done.\n')
nodeWeight = [];
NodeProbability = [];
fprintf('Finding critical nodes...')
for i = 1:size(spikeTrials{1},1)
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
fprintf('done.\n')
% Normalize node weights and find threshold for crit nodes
normNodeWeight = normalize(NodeWeight,'range');
normNodeProbability = normalize(NodeProbability,'range');
nodeWThresh = 0.6;
nodePThresh = 0.6;

% find crit nodes across each axis
critWNodes = find(normNodeWeight>=nodeWThresh);
critWNodeValue = normNodeWeight(normNodeWeight>=nodeWThresh);
critPNodes = find(normNodeProbability>=nodePThresh);
critPNodeValue = normNodeProbability(normNodeProbability>=nodePThresh);
% check membership for each 
critNodes = [];
if length(critWNodes)>length(critPNodes) %checks whick one is larger for membership assignment
    critNodes = critPNodes(ismember(critPNodes,critWNodes)); 
else 
    critNodes = critWNodes(ismember(critWNodes,critPNodes));
end
% Outputs
Ensemble.nodeWeight = NodeWeight;
Ensemble.nodeProbability = NodeProbability;
Ensemble.nodeWThresh = nodeWThresh;
Ensemble.nodePThresh = nodePThresh;
Ensemble.critNodes = critNodes;
