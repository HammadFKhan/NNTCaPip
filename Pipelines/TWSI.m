%% SI for temporal dynamics of TW ensembles
clear all
addpath(genpath('main'));
addpath(genpath('Pipeline'));

foldername = strcat(uigetdir(pwd,'Input Directory'),'\');
filetype = 'mat'; % type of files to be processed
% Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
if isempty(file),disp('No Mat files where detected in this directory!'), return; end % handing for incorrect files

numFile = length(file);
%%
for fileNum = 1:numFile
    filename = file(fileNum).name;
    load(filename)
    %     try
    %     %% Late Spike Corr Chunking
    %     trialType = trialData.responsiveTrials.lateSpikeTrials; %Assign trial id
    %     ensembleIdx = lateSpikeEnsemble.rankEnsembles(cellfun(@length,lateSpikeEnsemble.rankEnsembles)>2);
    %     ensembleIdx = vertcat(ensembleIdx{:});
    %     ensembleIdx = unique(ensembleIdx);
    %     corrChunkLate = chunkTrials(trialType,ensembleIdx,spikeTrials);
    %
    %     %% No Late Spike Corr Chunking
    %     trialType = trialData.responsiveTrials.noLateSpikeTrials; %Assign trial id
    %     ensembleIdx = nolateSpikeEnsemble.rankEnsembles(cellfun(@length,nolateSpikeEnsemble.rankEnsembles)>2);
    %     ensembleIdx = vertcat(ensembleIdx{:});
    %     ensembleIdx = unique(ensembleIdx);
    %
    %     corrChunknoLate = chunkTrials(trialType,ensembleIdx,spikeTrials);
    %     %% Do stats
    %     latecorrChunkmean = cellfun(@(x) mean(x,3),corrChunkLate,'UniformOutput',false);
    %     noLatecorrChunkmean = cellfun(@(x) mean(x,3),corrChunknoLate,'UniformOutput',false);
    %     % cell manipulation to 1-D
    %     lateSingleChunk = cell1D(latecorrChunkmean);
    %     noLateSingleChunk = cell1D(noLatecorrChunkmean);
    %     % Gather variables
    %     output.latecorrChunkmean = latecorrChunkmean;
    %     output.noLatecorrChunkmean = noLatecorrChunkmean;
    %     output.lateSingleChunk = lateSingleChunk;
    %     output.noLateSingleChunk = noLateSingleChunk;
    %     output.corrChunknoLate = corrChunknoLate;
    %     output.corrChunkLate = corrChunkLate;
    % Late vs no Late spike ensembles
    [lateSpikeEnsemble, nolateSpikeEnsemble] =...
        travelingWaveEnsemble(spikeTrials,trialData.responsiveTrials.lateSpikeTrials,trialData.responsiveTrials.noLateSpikeTrials,ROIcentroid,AverageImage);
    
    % manifold analysis and entropy
    lateSpikeEnsemble = ensembleMetric(lateSpikeEnsemble,AverageImage,ROIcentroid);
    nolateSpikeEnsemble = ensembleMetric(nolateSpikeEnsemble,AverageImage,ROIcentroid);
    
    % some statistics about these ensembles
    lateSpikeEnsemble = ensembleStat(lateSpikeEnsemble);
    nolateSpikeEnsemble = ensembleStat(nolateSpikeEnsemble);
    % Quantify Node Reactivation
    lateSpikeEnsemble = findCritNode(trialData.responsiveTrials.lateSpikeTrials,ROI,spikeTrials,lateSpikeEnsemble); % trials, ROI, spikeTrials, Spikes, Ensemble
    nolateSpikeEnsemble = findCritNode(trialData.responsiveTrials.noLateSpikeTrials,ROI,spikeTrials,nolateSpikeEnsemble);
    %% Plot some stuff
    %         for i = 1:11
    %             figure(1)
    %             subplot(4,3,i),imagesc(sum(corrChunkLate{i},3)),colormap(jet),caxis([0 1]),colorbar
    %             figure(2),
    %             NodeSize = 2;EdgeSize = 1;
    %             Connected_ROI = Connectivity_dice(latecorrChunkmean{i},.05);
    %             subplot(4,3,i),Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid(ensembleIdx,:),NodeSize,EdgeSize)
    %             figure(3)
    %             subplot(4,3,i),imagesc(sum(corrChunknoLate{i},3)),colormap(jet),caxis([0 1]),colorbar
    %             figure(4),
    %             NodeSize = 2;EdgeSize = 1;
    %             Connected_ROI = Connectivity_dice(noLatecorrChunkmean{i},.05);
    %             subplot(4,3,i),Cell_Map_Dice(AverageImage,Connected_ROI,ROIcentroid(ensembleIdx,:),NodeSize,EdgeSize)
    %         end
    %%
    figure(5),
    plot(smoothdata(output.lateSingleChunk)),hold on
    plot(smoothdata(output.noLateSingleChunk))
    legend('Late Spike', 'Weak Late Spike')
    
    %%
    if ~exist([file(fileNum).folder '\Analysis'],'dir')
        mkdir([file(fileNum).folder '\Analysis']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\Analysis'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\Analysis'],[file_name,'.mat']);
    save(savepath,'files', 'output');
    %     catch
    %         continue
    %     end
    clearvars -except file numFile fileNum filetype foldername
end
%% Local functions
% Late Spike Corr Chunking
function corrChunk = chunkTrials(trialType,ensembleIdx,spikeTrials)
chunkLength = 1:30:size(spikeTrials{1},2);% chunk the sequence into 1 sec intervals
%Analyze significant ensembles
try ensembleIdx = vertcat(ensembleIdx{:});ensembleIdx = unique(ensembleIdx);end
% Initialize
% corrChunk = cell(length(chunkLength),length(trialType));
corrChunk = cell(1,length(chunkLength)-1);
for trial = 1:length(trialType)
    disp(['Analyzing Trial: ' num2str(trial)])
    tempSpike = spikeTrials{trialType(trial)}; % Assign raster plot for specific trial
    for chunk = 1:length(chunkLength)-1
        corrChunk{trial}(:,:,chunk) = correlation_dice(tempSpike(ensembleIdx,chunkLength(chunk):chunkLength(chunk+1)-1));
    end
end
end

function singleChunk = cell1D(chunkMean)
singleChunk = cellfun(@(x) tril(x,-1), chunkMean, 'UniformOutput', false);
singleChunk = cell2mat(singleChunk);singleChunk = reshape(singleChunk,size(singleChunk,1),size(singleChunk,1),[]);
singleChunk(singleChunk==0) = NaN;
singleChunk = squeeze(nanmean(nanmean(singleChunk,1)));singleChunk(isnan(singleChunk)) = 0;
end
