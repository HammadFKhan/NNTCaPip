function centroidAnalysis(Ensemble)
%%
% Write centroid data
centroidData = [];
count = 1;
for i = 1:Ensemble.ensembleIndentified
    if ~isempty(Ensemble.ActivityCentroid{i})
        centroidData(count,1) = Ensemble.ActivityCentroid{i}(1);
        centroidData(count,2) = Ensemble.ActivityCentroid{i}(2);
        centroidData(count,3) = Ensemble.ActivityCentroidVariance{i}(1);
        centroidData(count,4) = Ensemble.ActivityCentroidVariance{i}(2);
        count = count+1;
    end
end

figure,errorbar(centroidData(:,1),centroidData(:,2),...
    centroidData(:,3),centroidData(:,3),centroidData(:,4),centroidData(:,4),'o'), hold on
legend
%% For loop for all Centroids

for i = 1:length(allCentroidData)
    figure(i),box off
    %normalize activity location (0,0)
    data = allCentroidData{i};
    data(:,1) = data(:,1)-min(data(:,1));
    data(:,2) = data(:,2)-min(data(:,2));
    errorbar(data(:,1),data(:,2),...
        data(:,3),data(:,3),data(:,4),data(:,4),'.'), hold on,box off
    axis([-200 250 -200 250])
end
%% Linear fitting everything
% mdl = fitlm(centroidFit(:,1),centroidFit(:,2));
% figure,plot(mdl,'marker','o','Color','r'),legend off,box off,title('Activity Centroid Area')
% mdl = fitlm(centroidFit(:,3),centroidFit(:,4));
% figure,plot(mdl,'marker','o','Color','b'),legend off,box off,title('Activity Centroid Variance')
% hyp1 = sqrt(centroidFit(:,1).^2 + centroidFit(:,2).^2);
% hyp2 = sqrt(centroidFit(:,3).^2 + centroidFit(:,4).^2);
% figure,plot(mdl,'marker','o','Color','g'),legend off,box off,title('Euclidean Centroid Activity vs. Variance')
for i = 1:length(allCentroidData)
%     hyp1 = sqrt(allCentroidData{i}(:,1).^2 + allCentroidData{i}(:,2).^2);
    
    hyp2 = sqrt(allCentroidData{i}(:,3).^2 + allCentroidData{i}(:,4).^2);
    mdl = fitlm(hyp1,hyp2);
    figure,plot(mdl,'marker','o','Color','k'),legend off,box off,hold on
    title(['R^2: ' num2str(mdl.Rsquared.Ordinary)])
end
mdl = fitlm(centroidFit(:,1),centroidFit(:,2))
figure,
h = plot(mdl);hold on;
delete(h(1))
for i = 1:5
    figure
    plot(allCentroidData{i}(:,3),allCentroidData{i}(:,4),'.'),axis([0 120 0 120]),box off
end
legend off
end