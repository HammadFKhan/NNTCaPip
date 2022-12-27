clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
%%
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S = load(fullfile(directory(i).folder,directory(i).name));
    StrName=fieldnames(S);
    StrName = StrName{1};
%     [P_c{i},I_c{i},H_c(i,1),H_e(i,1),sizeE{i},sizeEdge{i}] = EntropyParser(S.(StrName));
    [Ensemble(i).stat] = ensembleStat(S);
end 
%%
for i = 1:length(Ensemble)
   ensembleRecruitment(i) = Ensemble(i).stat.ensembleRecruitment;
   if isempty(Ensemble(i).stat.minRecruitment) 
       Ensemble(i).stat.minRecruitment = 0; 
   end
   minRecruitment(i) = Ensemble(i).stat.minRecruitment; 
   if isempty(Ensemble(i).stat.maxRecruitment)
       Ensemble(i).stat.maxRecruitment = 0; 
   end
   maxRecruitment(i) = Ensemble(i).stat.maxRecruitment; 
   if isempty(Ensemble(i).stat.ensembleSize)   
       Ensemble(i).stat.ensembleSize = 0; 
   end
   ensembleSizeMean(i) = mean(Ensemble(i).stat.ensembleSize);
end
%%
for i = 1:length(Ensemble)
    simI_tot{i} = mean(Ensemble(i).stat.sim_index);
    activityCentroidtot{i} = Ensemble(i).stat.activityCentroidVariance;
end
activityCentroidtot = vertcat(activityCentroidtot{:});
%%
data6 = test(25:end);
figure,boxplot(data,'PlotStyle','compact'),ylim([0 0.5]),box off,hold on
plot(1.1*ones(1,length(data)),data,'.');
%% Plot manifolds and centroids
% start = 1;stop = 5;
% start = 6;stop = 10;
% start = 11;stop = 16;
% start = 17;stop = 22;
% start = 23;stop = 29;
% start = 30;stop = 37;
%4 29 35
color = hsv((stop-start)+1);
m = (stop-start)+1;
hypVariance = [];
centroidData = [];
ensembleIntrinsicConnection = [];
count = 1;
figure
for ii = start:stop % Reference to recording session/animal
    manifold = rankedActivityCoords{ii}; % Eached rank coords is the distance between each cell pair that makes up an ensemble
    subplot(3,m,(1*m)-(m-count)),hold on,axis([0 600 0 600]),axis off
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
    
    subplot(3,m,(2*m)-(m-count)),hold on,axis([0 600 0 600]),axis off % Plots all the coordinates emsemble points embedded within the manifold
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
            x = manifold{i}(:,1);y = manifold{i}(:,2);
            scatter(x,y,4,color(count,:),'filled')
        end
    end
    
    % Centroid activity
    % Note that I forgot to resort this index to match the ranked activity
    % centroids! Hence, there is a mismatch between the centroid and
    % manifold/activity coordiantes. We need to recalculate here which is
    % taking the mean and std of the activity coordinates :)
    for i = 1:size(manifold,2)
        if size(manifold{i},1)>2
            %     centroidData = [activityCentroid{ii} activityCentroidVariance{ii}]; % Horzcat the centroid activity
            centroidData(i,:) = [mean(rankedActivityCoords{1,ii}{1,i}) std(rankedActivityCoords{1,ii}{1,i})];
            %     meanCentroidData = mean(centroidData);
            %     errorbar(meanCentroidData(:,1),meanCentroidData(:,2),...
            %         meanCentroidData(:,3),meanCentroidData(:,3),meanCentroidData(:,4),meanCentroidData(:,4),'.','Color',color(ii,:))
        end
    end
    numActivity = cellfun('size',rankedActivityCoords{1,ii},1); % Number of activity centroids
    numActivity(numActivity<3)= [];
    % Calculate adaptive threshold to display centroid data
    thresCum = cumsum(numActivity);
    thresh = find((thresCum./thresCum(end))>=0.8); % Plot centroid that encompases 80% of variability
    thresh = thresh(1);
    disp(thresh)
    subplot(3,m,(3*m)-(m-count)),hold on,axis([0 600 0 600]),axis off
    errorbar(centroidData(1:thresh,1),centroidData(1:thresh,2),...
        centroidData(1:thresh,3),centroidData(1:thresh,3),...
        centroidData(1:thresh,4),centroidData(1:thresh,4),'.','Color',color(count,:))
    hypVariance{count} = 2*sqrt(centroidData(:,3).^2+centroidData(:,4).^2); %Euclin distance across centroid variance for stats
    ensembleIntrinsicConnection{count} = numActivity'; % basically a measure of how many cells are
    centroidData = [];
    %recruited to the ensemble as a function of time. Should be similair to
    %ensemble recruitability
    count = count+1;
end

%% Line plot for variance
test = [];testStd = [];
test = cell2mat(cellfun(@mean,hypVariance,'Uniform',false));
testStd = cell2mat(cellfun(@std,hypVariance,'Uniform',false));
figure,hold on
boxplot(test,'plotStyle','compact')
plot(1.2*ones(length(test),1),test,'.')
figure,hold on
boxplot(testStd,'plotStyle','compact')
plot(1.2*ones(length(testStd),1),testStd,'.')
%% Centroid cross weeks
centroidActivityWeeks{6} = test;
centroidActivityVarianceWeeks{6} = testStd;
%% Plot ensemble intrinsic connection to overall variability
x = [];y = [];
figure,hold on,xlabel('Ensemble Intrinsic Connection'),ylabel('Centroid Variability')
x = vertcat(ensembleIntrinsicConnection{:});
y = vertcat(hypVariance{:});
plot(x,y,'k.'),axis([0 50 0 300])
figure,
mdl = fitlm(x,y);
plot(mdl),title(['R^2: ' num2str(mdl.Rsquared.Ordinary)]);axis([0 50 0 300])

%% Cross weeks
test = [];
test = cell2mat(cellfun(@mean,centroidActivityWeeks,'Uniform',false));
testStd = cell2mat(cellfun(@std,centroidActivityWeeks,'Uniform',false));
testSize = cellfun('size',centroidActivityWeeks,2);
t = testStd./sqrt(testSize);
errorbar(1:6,test,t),box off,xlim([ 0 7])
%% State Space Ratio (Variance/EnsembleRecruitment)
SSR{1} = y./x; %hypVariance/ensemble
%% Discritize sizes
xd = discretize(x,0:5:100); % 20 bins (points)
binEIC = xd*5;
for i = 1:20 % number of discritize bins
    binVariance{6,i} = y(xd==i); % map variance to binned size  
end
%% plot across weeks
figure,hold on
color = hsv(20);
for i = 1:20
    plot(test(:,i),'-o','Color',color(i,:))
end
%%
for i = 1:6
    figure,
    boxplot(SSR{i},'plotStyle','compact'),hold on
    plot(1.2*ones(length(SSR{i}),1),SSR{i},'.'),ylim([0 100]),box off
end
%%
color = hsv(6);
figure, hold on 
for i = 1:6
    plot(linspace((i-1)*10+1,i*10,length(SSR{i})),SSR{i},'.','Color',color(i,:))
end