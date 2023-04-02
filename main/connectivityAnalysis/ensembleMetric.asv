function Ensemble = ensembleMetric(Ensemble,AverageImage,ROIcentroid)
% ensembleVid(Ensemble,AverageImage,ROIcentroid,files);
% Displays ensemble overlay
[~,I] = sort(cellfun(@length,Ensemble.NodeList),'descend'); %sort max node size
rankEdges = Ensemble.NumEdges(:,I);
rankEnsembles = Ensemble.NodeList(:,I); 
[grad,~]=colorGradient([0 0 1] ,[0 0 0],5)
Ensemble.rankEnsembles = rankEnsembles;
if ~isempty(AverageImage) && ~isempty(ROIcentroid) && length(rankEnsembles)>2
    figure,
    for i = 1:5
        axis off
        color = jet(3);
        EnsembleMap(AverageImage,ROIcentroid,rankEnsembles{i},6,grad(i,:))
        set(gcf,'Position',[100 100 500 500])
        drawnow
        hold on
    end
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
Ensemble.informationEntropy = shannonEntropy(rankEnsembles,size(Ensemble.ensemble,1));
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
%% Create more plots to explain centroid variance/analysis

