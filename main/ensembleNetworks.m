function Ensemble = ensembleNetworks(Ensemble)
% Plot digraph with all active ensemble nodes
% Need to rerank Nodes
for ii = 1:Ensemble.ensembleIndentified
    diNodeList(:,2) = 1:length(Ensemble.NodeList{1,ii});
    diNodeList(:,1) = Ensemble.NodeList{1,ii};
    refList = Ensemble.Connected_ROI{1,ii}(:,1:2);
    for i = 1:size(diNodeList,1)
        refList(refList==diNodeList(i,1))= i;
    end
    G{ii} = graph(refList(:,1),refList(:,2));
    %%
    d = distances(G{ii});
    figure,p = plot(G{ii});
    bins = conncomp(G{ii});
    p.MarkerSize = 7;
%     p.EdgeColor = [0.5 0.5 .5];
    p.LineWidth = 1;
    p.NodeCData = bins;
    colormap(hsv(4))
    clear refList diNodeList
end
Ensemble.G = G;
