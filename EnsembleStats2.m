clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S = load(fullfile(directory(i).folder,directory(i).name));
    StrName=fieldnames(S);
    StrName = StrName{1};
%     [P_c{i},I_c{i},H_c(i,1),H_e(i,1),sizeE{i},sizeEdge{i}] = EntropyParser(S.(StrName));
    [ensembleSize{i},ensembleNum{i},ensembleRecruitment{i},minRecruitment{i},maxRecruitment{i},...
        rankedActivityCoords{i},activityCentroid{i},activityCentroidVariance{i}] = ensembleStat(S.(StrName));
end 
%%
figure,bar(cell2mat(ensembleRecruitment))
figure,bar(cell2mat(maxRecruitment))
figure,bar(cell2mat(minRecruitment))
figure,bar(cell2mat(ensembleNum))
for i = 1:33
    ensembleSizeMean(i) = mean(ensembleSize{i});
end
%%
data6 = test(25:end);
figure,boxplot(data,'PlotStyle','compact'),ylim([0 0.5]),box off,hold on
plot(1.1*ones(1,length(data)),data,'.');
%% Plot manifolds and centroids
for ii = 30:38
    xm = [];ym = [];
    manifold = rankedActivityCoords{ii};
    figure,hold on,axis([0 400 0 400])
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
            x = manifold{i}(:,1);y = manifold{i}(:,2);
            xm = [xm;x];ym = [ym;y];
            k = boundary(x,y);
            x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
            y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
            x2 = smoothdata(x1,'gaussian',50);
            y2 = smoothdata(y1,'gaussian',50);
            plot(x2,y2,'Color',[0.5 0.5 0.5 0.5])
        end
    end
    % Plot mean manifold for the boundry 
    x = xm;y = ym;
    k = boundary(x,y);
    x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
    y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
    x2 = smoothdata(x1,'gaussian',50);
    y2 = smoothdata(y1,'gaussian',50);
    plot(x2,y2,'Color',[0 0 0],'LineWidth',2)
    
    figure,hold on,axis([0 400 0 400])
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
            x = manifold{i}(:,1);y = manifold{i}(:,2);
            scatter(xm,ym,4,'k','filled')
        end
    end
end
%%
%%
figure,hold on,axis([0 400 0 400])
color = hsv(5);
for ii = 1:5
    centroidData = [activityCentroid{ii} activityCentroidVariance{ii}]; % Horzcat the centroid activity
    meanCentroidData = mean(centroidData);
    % errorbar(centroidData(:,1),centroidData(:,2),...
    %     centroidData(:,3),centroidData(:,3),centroidData(:,4),centroidData(:,4),'.')
    errorbar(meanCentroidData(:,1),meanCentroidData(:,2),...
        meanCentroidData(:,3),meanCentroidData(:,3),meanCentroidData(:,4),meanCentroidData(:,4),'.','Color',color(ii,:))
end