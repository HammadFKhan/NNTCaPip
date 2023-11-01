pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S(i).FOV = load(fullfile(directory(i).folder,directory(i).name));
end
%% Parse Ensemble stats
% Similarity
totalLsSim = [];
totalnLsSim = [];
for i = 1:L
    LsSim = S(i).FOV.lateSpikeEnsemble.LsSim;
    nLsSim = S(i).FOV.nolateSpikeEnsemble.nLsSim;
    totalLsSim{i} = mean(LsSim)';
    totalnLsSim{i} = mean(nLsSim)';
end
totalLsSim = vertcat(totalLsSim{:});
totalnLsSim = vertcat(totalnLsSim{:});
figure,boxplot(totalLsSim),box off,ylim([0 0.8]),title('Ls'),set(gca,'TickDir','out');
figure,boxplot(totalnLsSim),box off,ylim([0 0.8]),title('nLs'),set(gca,'TickDir','out');


%% Neurons per ensemble
totalLsEnsembleSize = [];
totalnLsEnsembleSize = [];

for i = 1:L
    LsEnsembleSize = S(i).FOV.lateSpikeEnsemble.rankEnsembles;
    LsEnsembleSize = cellfun('size',LsEnsembleSize,1);
    nLsEnsembleSize = S(i).FOV.nolateSpikeEnsemble.rankEnsembles;
    nLsEnsembleSize = cellfun('size',nLsEnsembleSize,1);
    totalLsEnsembleSize{i} = LsEnsembleSize(LsEnsembleSize>2)';
    totalnLsEnsembleSize{i} = nLsEnsembleSize(nLsEnsembleSize>2)';
end
totalLsEnsembleSize = vertcat(totalLsEnsembleSize{:});
totalnLsEnsembleSize = vertcat(totalnLsEnsembleSize{:});
figure,boxplot(totalLsEnsembleSize),box off,ylim([0 30]),title('Ls'),set(gca,'TickDir','out');
figure,boxplot(totalnLsEnsembleSize),box off,ylim([0 30]),title('nLs'),set(gca,'TickDir','out');

%% Connections per ensemble
totalLsConnections = [];
totalnLsConnections = [];

for i = 1:L
    LsConnections = S(i).FOV.lateSpikeEnsemble.rankEdges;
    LsConnections = vertcat(LsConnections{:});
    nLsConnections = S(i).FOV.nolateSpikeEnsemble.rankEdges;
    nLsConnections = vertcat(nLsConnections{:});
    totalLsConnections{i} = mean(LsConnections(LsConnections>0));
    totalnLsConnections{i} = mean(nLsConnections(nLsConnections>5));
end
totalLsConnections = vertcat(totalLsConnections{:});
totalnLsConnections = vertcat(totalnLsConnections{:});
figure,boxplot(totalLsConnections),box off,ylim([0 30]),title('Ls'),set(gca,'TickDir','out');
figure,boxplot(totalnLsConnections),box off,ylim([0 30]),title('nLS'),set(gca,'TickDir','out');

%% Ensemble Recruitment
for i = 1:L
    lsRecruitment(i) = S(i).FOV.lateSpikeEnsemble.ensembleRecruitment;
    nlsRecruitment(i) = S(i).FOV.nolateSpikeEnsemble.ensembleRecruitment;
end
figure,boxplot(lsRecruitment),box off, ylim([0 0.1]),title('Ls'),set(gca,'TickDir','out');
figure,boxplot(nlsRecruitment),box off, ylim([0 0.1]),title('nLs'),set(gca,'TickDir','out');

%% Ensemble Weight
for i = 1:L
    lsrankEnsembles = S(i).FOV.lateSpikeEnsemble.rankEnsembles;
    sizeE = cellfun(@size,lsrankEnsembles,'UniformOutput',false);
    sizeE = cell2mat(sizeE);
    sizeE(sizeE==1) = [];
    sizeE = sizeE/max(sizeE);
    lsSizeE{i} = sizeE;
    
    nlsrankEnsembles = S(i).FOV.nolateSpikeEnsemble.rankEnsembles;
    sizeE = cellfun(@size,nlsrankEnsembles,'UniformOutput',false);
    sizeE = cell2mat(sizeE);
    sizeE(sizeE==1) = [];
    sizeE = sizeE/max(sizeE);
    nlsSizeE{i} = sizeE;
    figure
    plot(lsSizeE{i}),hold on, box off
    plot(nlsSizeE{i})
end
t = sort(horzcat(lsSizeE{:}),'descend');
t1 = sort(horzcat(nlsSizeE{:}),'descend');
figure,plot(t),box off,axis([0 200 .3 1]),hold on
plot(t1)

%% Neuron IDs for each FOV ensemble
for i = 1:L
    FOV(i).name = directory(i).name;
    FOV(i).lateSpike.neuronID = vertcat(S(i).FOV.lateSpikeEnsemble.rankEnsembles{1:3});
    FOV(i).nolateSpike.neuronID = vertcat(S(i).FOV.nolateSpikeEnsemble.rankEnsembles{1:3});
end


%%
totalActivityCentroid = horzcat(totalActivityCentroid,S(4).FOV.ActivityCentroid);
totalConnectedROI = horzcat(totalConnectedROI,S(4).FOV.Connected_ROI);
totalROIcentroid = vertcat(totalROIcentroid, S(4).FOV.ROIcentroid);
totalSpikes = horzcat(S(1).FOV.Spikes,S(2).FOV.Spikes);
%%
figure,
NodeSize = 1;EdgeSize = 1;
for i = 1:length(lateSpikeTrials(lateSpikeTrials<123))
    subplot(9,9,i),Cell_Map_Dice(AverageImage,totalConnectedROI{lateSpikeTrials(i)},totalROIcentroid,NodeSize,EdgeSize)
end
suptitle('Late Spike')

figure,
NodeSize = 1;EdgeSize = 1;
for i = 1:length(nolateSpikeTrials(nolateSpikeTrials<123))
    subplot(6,6,i),Cell_Map_Dice(AverageImage,totalConnectedROI{nolateSpikeTrials(i)},totalROIcentroid,NodeSize,EdgeSize)
end
suptitle('No Late Spike')


