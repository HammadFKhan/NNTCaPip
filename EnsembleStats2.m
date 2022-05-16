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
color = hsv(5);
for ii = 1:5 % Reference to recording session/animal
    manifold = rankedActivityCoords{ii}; % Eached rank coords is the distance between each cell pair that makes up an ensemble
    figure,hold on,axis([0 400 0 400])
    for i = 1:size(manifold,2) % Basically how many ensembles there are
        if size(manifold{i},1)>2
            x = manifold{i}(:,1);y = manifold{i}(:,2);
            k = boundary(x,y);
            x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
            y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
            x2 = smoothdata(x1,'gaussian',50);
            y2 = smoothdata(y1,'gaussian',50);
            plot(x2,y2,'Color',[0.5 0.5 0.5 0.5])
        end
    end
    % Plot mean manifold for the boundry 
%     k = boundary(x,y);
%     x1 = interp1(1:length(x(k)),x(k),1:0.05:length(x(k)),'pchip');
%     y1 = interp1(1:length(y(k)),y(k),1:0.05:length(y(k)),'pchip');
%     x2 = smoothdata(x1,'gaussian',50);
%     y2 = smoothdata(y1,'gaussian',50);
%     plot(x2,y2,'Color',[0 0 0],'LineWidth',2)
    
    figure,hold on,axis([0 400 0 400]) % Plots all the coordinates emsemble points embedded within the manifold
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
            x = manifold{i}(:,1);y = manifold{i}(:,2);
            scatter(x,y,4,color(ii,:),'filled')
        end
    end
    
    % Centroid activity
    % Note that I forgot to resort this index to match the ranked activity
    % centroids! Hence, there is a mismatch between the centroid and
    % manifold/activity coordiantes. We need to recalculate here which is
    % taking the mean and std of the activity coordinates :)
    figure,hold on,axis([0 400 0 400])
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
        %     centroidData = [activityCentroid{ii} activityCentroidVariance{ii}]; % Horzcat the centroid activity
        
        % adaptively find threshold of centoriod activity that encompases
        % 80% of ensemble variability
        
        centroidData(i,:) = [mean(rankedActivityCoords{1,ii}{1,i}) std(rankedActivityCoords{1,ii}{1,i})];
        %     meanCentroidData = mean(centroidData);
        %     errorbar(meanCentroidData(:,1),meanCentroidData(:,2),...
        %         meanCentroidData(:,3),meanCentroidData(:,3),meanCentroidData(:,4),meanCentroidData(:,4),'.','Color',color(ii,:))
        end
    end
    errorbar(centroidData(:,1),centroidData(:,2),...
            centroidData(:,3),centroidData(:,3),centroidData(:,4),centroidData(:,4),'.','Color',color(ii,:))
    hypVariance{ii} = 2*sqrt(centroidData(:,3).^2+centroidData(:,4).^2); %Euclin distance across centroid variance for stats
    centroidData = [];
end

%% Line plot for variance
test = cell2mat(cellfun(@mean,hypVariance,'Uniform',false));
testStd = cell2mat(cellfun(@std,hypVariance,'Uniform',false));



%% Centroid activity
% Note that I forgot to resort this index to match the ranked activity
% centroids! Hence, there is a mismatch between the centroid and
% manifold/activity coordiantes. We need to recalculate here which is
% taking the mean and std of the activity coordinates :)
figure,hold on,axis([0 400 0 400])
color = hsv(5);
for ii = 5
%     centroidData = [activityCentroid{ii} activityCentroidVariance{ii}]; % Horzcat the centroid activity
centroidData = [mean(rankedActivityCoords{1,ii}{1,ii}) std(rankedActivityCoords{1,ii}{1,ii})];
%     meanCentroidData = mean(centroidData);
    errorbar(centroidData(:,1),centroidData(:,2),...
        centroidData(:,3),centroidData(:,3),centroidData(:,4),centroidData(:,4),'.')
%     errorbar(meanCentroidData(:,1),meanCentroidData(:,2),...
%         meanCentroidData(:,3),meanCentroidData(:,3),meanCentroidData(:,4),meanCentroidData(:,4),'.','Color',color(ii,:))
end