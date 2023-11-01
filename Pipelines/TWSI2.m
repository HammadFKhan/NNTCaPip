set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'));
std_threshold = 6;
static_threshold = 0;
Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);


spikeTrials = [];
trialLength = 360;
for i = 1:size(Spikes,2)/trialLength
    spikeTrials{i} = Spikes(:,((i-1)*trialLength+1):i*trialLength);
    DeltaTrials(:,:,i) = DeltaFoverF(:,((i-1)*trialLength+1):i*trialLength);
end

% check for edge case when calcium and trial numbers don't match
if length(spikeTrials)<trialData.responsiveTrials.trialNum(end)
    trialData.responsiveTrials.lateSpikeTrialsOld = trialData.responsiveTrials.lateSpikeTrials;
    trialData.responsiveTrials.lateSpikeTrials = trialData.responsiveTrials.lateSpikeTrials...
        (trialData.responsiveTrials.lateSpikeTrials<=length(spikeTrials));
    
    trialData.responsiveTrials.noLateSpikeTrialsOld = trialData.responsiveTrials.noLateSpikeTrials;
    trialData.responsiveTrials.noLateSpikeTrials = trialData.responsiveTrials.noLateSpikeTrials...
        (trialData.responsiveTrials.noLateSpikeTrials<=length(spikeTrials));
end


% Late vs no Late spike ensembles
[lateSpikeEnsemble, nolateSpikeEnsemble] =...
    travelingWaveEnsemble(spikeTrials,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROIcentroid,AverageImage);

% manifold analysis and entropy
lateSpikeEnsemble = ensembleMetric(lateSpikeEnsemble,AverageImage,ROIcentroid);
nolateSpikeEnsemble = ensembleMetric(nolateSpikeEnsemble,AverageImage,ROIcentroid);

% some statistics about these ensembles
lateSpikeEnsemble = ensembleStat(lateSpikeEnsemble);
nolateSpikeEnsemble = ensembleStat(nolateSpikeEnsemble);
close all
%%
clear all
addpath(genpath('main'));
addpath(genpath('Pipeline'));

foldername = strcat(uigetdir(pwd,'Input Directory'),'\');
filetype = 'mat'; % type of files to be processed
% Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
if isempty(file),disp('No Mat files where detected in this directory!'), return; end % handing for incorrect files

numFile = length(file);
for fileNum = 1:numFile
    filename = file(fileNum).name;
    load(filename)
    % Now do the same thing but with truncated sensory driven trials
    spikeTrialsTemp = cellfun(@(x) x(:,70:179),spikeTrials,'UniformOutput',false);
    [lateSpikeEnsemblesensory, nolateSpikeEnsemblesensory] =...
        travelingWaveEnsemble(spikeTrialsTemp,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROIcentroid,AverageImage);
    
    % manifold analysis and entropy
    lateSpikeEnsemblesensory = ensembleMetric(lateSpikeEnsemblesensory,AverageImage,ROIcentroid);
    nolateSpikeEnsemblesensory = ensembleMetric(nolateSpikeEnsemblesensory,AverageImage,ROIcentroid);
    
    % some statistics about these ensembles
    lateSpikeEnsemblesensory = ensembleStat(lateSpikeEnsemblesensory);
    nolateSpikeEnsemblesensory = ensembleStat(nolateSpikeEnsemblesensory);
    close all
    % Now do the same thing but with truncated pre-sensory driven trials
    spikeTrialsTemp = cellfun(@(x) x(:,1:60),spikeTrials,'UniformOutput',false);
    [lateSpikeEnsemblepre, nolateSpikeEnsemblepre] =...
        travelingWaveEnsemble(spikeTrialsTemp,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROIcentroid,AverageImage);
    
    % manifold analysis and entropy
    lateSpikeEnsemblepre = ensembleMetric(lateSpikeEnsemblepre,AverageImage,ROIcentroid);
    nolateSpikeEnsemblepre = ensembleMetric(nolateSpikeEnsemblepre,AverageImage,ROIcentroid);
    
    % some statistics about these ensembles
    lateSpikeEnsemblepre = ensembleStat(lateSpikeEnsemblepre);
    nolateSpikeEnsemblepre= ensembleStat(nolateSpikeEnsemblepre);
    close all
    
    simM = [];
    Connected_ROI = [];
    senLs = [];
    sim_indexsen = [];
    sim_indexpre = [];
    
    [~,sim_indexsen] = cosine_similarity(lateSpikeEnsemblesensory.Spikes,10);
    [~,sim_indexpre] = cosine_similarity(lateSpikeEnsemblepre.Spikes,10);
    senLs = mean(sim_indexsen(sim_indexsen>0.1),'all'); %mean center
    preLs = mean(sim_indexpre,'all'); %mean center
    lateSpikeEnsemblesensory.LsSim = senLs;
    lateSpikeEnsemblepre.LsSim = preLs;
    
    
    Connected_ROI = {};
    simM = [];
    simValue = [];
    nLs = [];
    sim_indexsen = [];
    sim_indexpre = [];
    % % figure,
    [~,sim_indexsen] = cosine_similarity(nolateSpikeEnsemblesensory.Spikes,10);
    [~,sim_indexpre] = cosine_similarity(nolateSpikeEnsemblepre.Spikes,10);
    sennLs = mean(sim_indexsen(sim_indexsen>0),'all');
    prenLs = mean(sim_indexpre,'all'); %mean center
    
    nolateSpikeEnsemblesensory.nLsSim = sennLs;
    nolateSpikeEnsemblepre.nLsSim = prenLs;
    %%
    if ~exist([file(fileNum).folder '\EnsembleTimeAnalysis'],'dir')
        mkdir([file(fileNum).folder '\EnsembleTimeAnalysis']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\EnsembleTimeAnalysis'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\EnsembleTimeAnalysis'],[file_name,'.mat']);
    save(savepath,'files','ROIcentroid','lateSpikeEnsemble','lateSpikeEnsemblesensory','lateSpikeEnsemblepre',...
        'nolateSpikeEnsemble','nolateSpikeEnsemblesensory','nolateSpikeEnsemblepre');
    %     catch
    %         continue
    %     end
    clearvars -except file numFile fileNum filetype foldername
end

%% Plot stuff
figure,
subplot(1,3,1),imagesc(lateSpikeEnsemble.sim_index),colormap(hot)
subplot(1,3,2),imagesc(lateSpikeEnsemblesensory.sim_index),colormap(hot)
subplot(1,3,3),imagesc(lateSpikeEnsemblepre.sim_index),colormap(hot)
set(gcf,'Position',[100 100 1500 350])

figure,
subplot(1,3,1),imagesc(nolateSpikeEnsemble.sim_index),colormap(hot)
subplot(1,3,2),imagesc(nolateSpikeEnsemblesensory.sim_index),colormap(hot)
subplot(1,3,3),imagesc(nolateSpikeEnsemblepre.sim_index),colormap(hot)
set(gcf,'Position',[100 100 1500 350])

t1 = lateSpikeEnsemble.sim_index; t1(t1==0) = NaN;t1 = nanmean(t1)';
t2 = lateSpikeEnsemblesensory.sim_index; t2(t2==0) = NaN;t2 = nanmean(t2)';
t3 = lateSpikeEnsemblepre.sim_index; t3(t3==0) = NaN;t3 = nanmean(t3)';

t4 = nolateSpikeEnsemble.sim_index; t4(t4==0) = NaN;t4 = nanmean(t4)';
t5 = nolateSpikeEnsemblesensory.sim_index; t5(t5==0) = NaN;t5 = nanmean(t5)';
t6 = nolateSpikeEnsemblepre.sim_index; t6(t6==0) = NaN;t6 = nanmean(t6)';

figure,
for i = 2:5
    axis off
    color = jet(3);
    EnsembleMap(AverageImage,ROIcentroid,nolateSpikeEnsemblepre.rankEnsembles{i},5,[0 0 1])
    set(gcf,'Position',[100 100 500 500])
    drawnow
    hold on
end
%-------------------------------------------------------%
%-------------------------------------------------------%
%-------------------------------------------------------%
%% Stats across animals
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S(i).FOV = load(fullfile(directory(i).folder,directory(i).name));
end
%% Parse Ensemble stats
% Similarity index
totalLsSim = [];
totalLsSimSen = [];
totalLsSimPre = [];
totalnLsSim = [];
totalnLsSimSen = [];
totalnLsSimPre = [];
for i = 1:L
    %     LsSim = S(i).FOV.lateSpikeEnsemble.sim_index;
    [~,LsSimSen] = cosine_similarity(S(i).FOV.lateSpikeEnsemblesensory.Spikes,10);
    [~,LsSimPre] = cosine_similarity(S(i).FOV.lateSpikeEnsemblepre.Spikes,10);
    
    [~,nLsSimSen] = cosine_similarity(S(i).FOV.nolateSpikeEnsemblesensory.Spikes,10);
    [~,nLsSimPre] = cosine_similarity(S(i).FOV.nolateSpikeEnsemblepre.Spikes,10);
    
    %     subplot(5,4,count),imagesc(sim_index),colormap(jet)
%         LsSimSen = S(i).FOV.lateSpikeEnsemblesensory.LsSim;
%         LsSimPre = S(i).FOV.lateSpikeEnsemblepre.LsSim;
% %         nLsSim = S(i).FOV.nolateSpikeEnsemble.nLsSim;
%         nLsSimSen = S(i).FOV.nolateSpikeEnsemblesensory.nLsSim;
%         nLsSimPre = S(i).FOV.nolateSpikeEnsemblepre.nLsSim;
%         totalLsSim{i} = mean(LsSim(LsSim>0))';
%         totalnLsSim{i} = mean(nLsSim(nLsSim>0))';
    totalLsSimSen{i} = mean(LsSimSen(LsSimSen>0.15))';
    totalnLsSimSen{i} = mean(nLsSimSen(nLsSimSen>0))';
    totalLsSimPre{i} = mean(LsSimPre)';
    totalnLsSimPre{i} = mean(nLsSimPre)';
end
%%
totalLsSim = vertcat(totalLsSim{:});
totalLsSimSen = vertcat(totalLsSimSen{:});
totalLsSimPre = vertcat(totalLsSimPre{:});
totalnLsSim = vertcat(totalnLsSim{:});
totalnLsSimSen = vertcat(totalnLsSimSen{:});
totalnLsSimPre = vertcat(totalnLsSimPre{:});

figure,customBoxplot([totalLsSim,totalnLsSim,totalLsSimSen,totalnLsSimSen,totalLsSimPre,totalnLsSimPre])
%%
ensembleAll = vertcat(S(i).FOV.lateSpikeEnsemblepre.rankEnsembles{:});
ensembleAll(ensembleAll==0) = [];
ensembleAll = unique(ensembleAll);
corr = correlation_dice(S(i).FOV.lateSpikeEnsemblepre.ensemble(ensembleAll,:));
figure,imagesc(corr),colormap(hot)
ensembleAll = vertcat(S(i).FOV.nolateSpikeEnsemblepre.rankEnsembles{:});
ensembleAll(ensembleAll==0) = [];
ensembleAll = unique(ensembleAll);
corr = correlation_dice(S(i).FOV.nolateSpikeEnsemblepre.ensemble(ensembleAll,:));
figure,imagesc(corr),colormap(hot)
%% Neurons per ensemble
totalLsEnsembleSize = [];
totalnLsEnsembleSize = [];

for i = 1:L
    LsEnsembleSize = S(i).FOV.lateSpikeEnsemblesensory.rankEnsembles;
    LsEnsembleSize = cellfun('size',LsEnsembleSize,1);
    nLsEnsembleSize = S(i).FOV.nolateSpikeEnsemblesensory.rankEnsembles;
    nLsEnsembleSize = cellfun('size',nLsEnsembleSize,1);
    LsEnsembleSizepre = S(i).FOV.lateSpikeEnsemblepre.rankEnsembles;
    LsEnsembleSizepre = cellfun('size',LsEnsembleSizepre,1);
    nLsEnsembleSizepre = S(i).FOV.nolateSpikeEnsemblepre.rankEnsembles;
    nLsEnsembleSizepre = cellfun('size',nLsEnsembleSizepre,1);
    totalLsEnsembleSize{i} = LsEnsembleSize(LsEnsembleSize>2)';
    totalnLsEnsembleSize{i} = nLsEnsembleSize(nLsEnsembleSize>2)';
    totalLsEnsembleSizepre{i} = LsEnsembleSizepre(LsEnsembleSizepre>2)';
    totalnLsEnsembleSizepre{i} = nLsEnsembleSizepre(nLsEnsembleSizepre>2)';
end

totalLsEnsembleSize = vertcat(totalLsEnsembleSize{:});
totalnLsEnsembleSize = vertcat(totalnLsEnsembleSize{:});
totalLsEnsembleSizepre = vertcat(totalLsEnsembleSizepre{:});
totalnLsEnsembleSizepre = vertcat(totalnLsEnsembleSizepre{:});

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
%% PCA instability from trial data
Magnitude = cell(32,1);
for i = 1:32
    win = size(lateSpikeEnsemble.Spikes,2)-(i-1)*100;
    Spikes = lateSpikeEnsemble.Spikes(:,1:win);
    Magnitude{i} = PCAvariability(Spikes);
end
%%
t = cellfun(@(x) nanmean(x,1),Magnitude1,'UniformOutput',false);










%%
function Magnitude = PCAvariability(Spikes)
disp(['Shuffling data ' num2str(1000) ' times to find optimal ensemble cutoff']);
shufSpikes = tempShuffle(Spikes,1000);
[coactive_cells,~] = coactive_index(Spikes,length(Spikes)/5);
[shufcoactive_cells,~] = coactive_index(shufSpikes,length(Spikes)/5);
bin = ceil(max(coactive_cells)*100);
figure,hold on,bar(coactive_cells),bar(shufcoactive_cells);
% figure,hold on,histogram(coactive_cells,100),histogram(shufcoactive_cells,100);
ensembleWin = 5;
ensembleCan = ensembleWin*find(coactive_cells>(2.5*std(shufcoactive_cells)+mean(shufcoactive_cells))); % 99% distribution threshold
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
    [vectorized,sim_index] = cosine_similarity(ensemble(:,1:factorCorrection),ensembleWin);
end
fprintf('done\n')
avgSim = mean(mean(sim_index,2));
X = sim_index - avgSim*ones(1,size(sim_index,2));
[U,S,V] = svd(X,'econ');
plotind = 1;
try
    figure,
    for r = [1:25]
        Xapprox = U(:,r)*S(r,r)*V(:,r)';
        subplot(5,5,plotind), plotind = plotind+1;
        imagesc(Xapprox),colormap(flip(gray)),caxis([0 0.15]);
    end
catch
end
% Weighted Ensemble
ensembleSimilarity = zeros(1,size(S,1));
for i = 1:size(S,1)
    Xapprox = U(:,i)*S(i,i)*V(:,i)';
    ensembleSimilarity(i) = mean(Xapprox(Xapprox>0))+avgSim;
end
% Singular Values
% figure,subplot(1,3,1)
% semilogx(ensembleSimilarity,'k','Linewidth',2),grid on;hold on;
% xlabel('r'); ylabel('Ensemble Similarity')
% set(gca,'FontSize',14);
% subplot(1,3,2)
% semilogx(diag(S),'k','Linewidth',2),grid on;hold on;
% xlabel('r'); ylabel('Singular Value, \Sigma_r')
% set(gca,'FontSize',14);
% subplot(1,3,3)
% plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2),grid on;
% xlabel('r'); ylabel('Cumulative Energy');
% set(gca,'FontSize',14);
Magnitude(1) = nanmean(ensembleSimilarity);
Magnitude(2) = nanstd(ensembleSimilarity);

%-----------Now use the shuffled version-----------------%
Spikes = shufSpikes;
ensembleWin = 5;
ensembleCan = ensembleWin*find(shufcoactive_cells>(2.5*std(shufcoactive_cells)+mean(shufcoactive_cells))); % 99% distribution threshold
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
    [vectorized,sim_index] = cosine_similarity(ensemble(:,1:factorCorrection),ensembleWin);
end
fprintf('done\n')
avgSim = mean(mean(sim_index,2));
X = sim_index - avgSim*ones(1,size(sim_index,2));
[U,S,V] = svd(X,'econ');
plotind = 1;
try
    figure,
    for r = [1:25]
        Xapprox = U(:,r)*S(r,r)*V(:,r)';
        subplot(5,5,plotind), plotind = plotind+1;
        imagesc(Xapprox),colormap(flip(gray)),caxis([0 0.15]);
    end
catch
end
% Weighted Ensemble
ensembleSimilarity = zeros(1,size(S,1));
for i = 1:size(S,1)
    Xapprox = U(:,i)*S(i,i)*V(:,i)';
    ensembleSimilarity(i) = mean(Xapprox(Xapprox>0))+avgSim;
end
% Singular Values
% figure,subplot(1,3,1)
% semilogx(ensembleSimilarity,'k','Linewidth',2),grid on;hold on;
% xlabel('r'); ylabel('Ensemble Similarity')
% set(gca,'FontSize',14);
% subplot(1,3,2)
% semilogx(diag(S),'k','Linewidth',2),grid on;hold on;
% xlabel('r'); ylabel('Singular Value, \Sigma_r')
% set(gca,'FontSize',14);
% subplot(1,3,3)
% plot(cumsum(diag(S))/sum(diag(S)),'k','LineWidth',2),grid on;
% xlabel('r'); ylabel('Cumulative Energy');
% set(gca,'FontSize',14);
Magnitude(3) = nanmean(ensembleSimilarity);
Magnitude(4) = nanstd(ensembleSimilarity);
end



