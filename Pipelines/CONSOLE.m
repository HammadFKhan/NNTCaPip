%% NNT Calcium Pipeline
% Github Version 4.1
ds_filename = intanPreprocessingImaging;
data = matfile(ds_filename);
fpath = data.fpath;
savepath = fullfile(fpath,['loadme','.mat']);
save(savepath,'ds_filename');
clearvars -except ds_filename
%%
data = matfile(ds_filename);
parameters.experiment = 'cue'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.windowBeforePull = 0.5; % in seconds
parameters.windowAfterPull = 4.5; % in seconds
parameters.windowBeforeCue = 0.5; % in seconds
parameters.windowAfterCue = 4.5; % in seconds
parameters.windowBeforeMI = 0.5; % in seconds
parameters.windowAfterMI = 4.5; % in seconds
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
amplifierTime = downsample(data.amplifierTime,round(5000/parameters.Fs),1); % time in seconds
[Behaviour] = readLever(parameters,amplifierTime);
[IntanBehaviour] = readLeverIntanImaging(parameters,amplifierTime,data.analogChannels(1,:),data.digitalChannels,Behaviour,1);
IntanBehaviour.parameters = parameters;
% Calculate ITI time for trials and reward/no reward sequence
temp1 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueHitTrace);
temp1 = vertcat(temp1,ones(1,IntanBehaviour.nCueHit)); %  write 1 for reward given
temp2 = arrayfun(@(x) x.LFPtime(1), IntanBehaviour.cueMissTrace);
temp2 = vertcat(temp2,zeros(1,IntanBehaviour.nCueMiss)); %  write 0 for no reward given
temp = [temp1,temp2];
[~,idx] = sort(temp(1,:)); %sort by occurance
IntanBehaviour.ITI = temp(:,idx);
%% Remove ROIs
if exist('badComponents','var') && ~exist('badComFlag','var')
    [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
        removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(156));
    badComFlag = 1;
end
% Fix centroids
ROIcentroid = [];
for i = 1:length(ROI)
    blah = vertcat(ROI{i}{:});
    ROIcentroid(i,:) = floor(mean(blah,1));
end
%% preAnalysis (setup paths and grab sutter metadata)
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'));
[fileInfo] = CaImAnReadIn([]); %grabs metadata
%% Volume Detection (volumetric concatenation
if volumeCa
    l5Idx = 1:size(DeltaFoverF,1);
    l5DeltaFoverF = DeltaFoverF;
end
%% Spike detection from dF/F
std_threshold = 3;      % from Carrilo-Reid and Jordan Hamm's papers
static_threshold = .01;
Spikes = rasterizeDFoF(DeltaFoverF,std_threshold,static_threshold);
figure('Name','Spiking Raster');Show_Spikes(Spikes);
%%
%Excude inactive cells
% numSpikes = sum(Spikes,2);
% keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
% Spikes = Spikes(keepSpikes,:);
[coactive_cells,detected_spikes] = coactive_index(Spikes,size(Spikes,2));
cell_count = length(ROI);
caFR = 15.0282; %% 30.0647;
time = time_adjust(size(DeltaFoverF,2),caFR);
for i = 1:size(DeltaFoverF,1)
    calcium_avg{i} = STA(DeltaFoverF(i,:),2,250);%std, window (frames)
end
% Pairwise Velocity Analysis
% velocityPairwise(VR_data,Spikes)
%% Behavioural analysis of calcium and spikes
Calcium = leverCaModulation(l23DeltaFoverF,l23Spikes,IntanBehaviour,time);
%% Ensemble Analysis
idx = size(Calcium.hit.Spikes);
hitSpikes = reshape(Calcium.hit.Spikes,idx(1),idx(2)*idx(3));
idx = size(Calcium.miss.Spikes);
missSpikes = reshape(Calcium.miss.Spikes,idx(1),idx(2)*idx(3));
idx = size(Calcium.MIFA.Spikes);
FASpikes = reshape(Calcium.MIFA.Spikes,idx(1),idx(2)*idx(3));
% Hit Ensemble
factorCorrection = 5*floor(size(hitSpikes,2)/5); % Correct for frame size aquisition
hitEnsemble = ensembleAnalysis2(hitSpikes(:,1:factorCorrection),ROIcentroid);
hitEnsemble = ensembleMetric(hitEnsemble,AverageImage,ROIcentroid);
hitEnsemble = ensembleStat(hitEnsemble);
% Miss Ensemble
factorCorrection = 5*floor(size(missSpikes,2)/5); % Correct for frame size aquisition
missEnsemble = ensembleAnalysis2(missSpikes(:,1:factorCorrection),ROIcentroid);
missEnsemble = ensembleMetric(missEnsemble,AverageImage,ROIcentroid);
missEnsemble = ensembleStat(missEnsemble);
% FA Ensemble
factorCorrection = 5*floor(size(FASpikes,2)/5); % Correct for frame size aquisition
FAEnsemble = ensembleAnalysis2(FASpikes(:,1:factorCorrection),ROIcentroid);
FAEnsemble = ensembleMetric(FAEnsemble,AverageImage,ROIcentroid);
FAEnsemble = ensembleStat(FAEnsemble);
%%
timeWin = [size(hitEnsemble.ensemble,2),size(missEnsemble.ensemble,2),size(FAEnsemble.ensemble,2)];
timeWin = 1:floor(min(timeWin)/10)*10;
tEnsemble = horzcat(hitEnsemble.ensemble(:,timeWin),missEnsemble.ensemble(:,timeWin),FAEnsemble.ensemble(:,timeWin));
[tVectorized,tsimIndex] = cosine_similarity(tEnsemble,5);
[X] = featureProject(tsimIndex, 1, 0);
%%
hitWin = 1:size(tsimIndex,2)/3;
missWin = 1+hitWin(end):(hitWin(end)+size(tsimIndex,2)/3);
FAWin = 1+missWin(end):(missWin(end)+size(tsimIndex,2)/3);
figure,
scatter3(X(hitWin,1),X(hitWin,2),X(hitWin,3),20,'filled'),hold on
scatter3(X(missWin,1),X(missWin,2),X(missWin,3),20,'filled','r')
scatter3(X(FAWin,1),X(FAWin,2),X(FAWin,3),20,'filled','k')
%%
figure,bar([hitEnsemble.ensembleNum,missEnsemble.ensembleNum,FAEnsemble.ensembleNum]),box off
set(gca,'TickDir','out','fontsize',14),ylabel('# of Ensembles')
figure,bar([mean(hitEnsemble.ensembleSize),mean(missEnsemble.ensembleSize),mean(FAEnsemble.ensembleSize)]),box off
set(gca,'TickDir','out','fontsize',14),ylabel('# Neurons/Ensembles')

%%
hitWin = 1:size(tsimIndex,2)/3;
missWin = 1+hitWin(end):(hitWin(end)+size(tsimIndex,2)/3);
FAWin = 1+missWin(end):(missWin(end)+size(tsimIndex,2)/3);
ensembleStab = [mean(tsimIndex(hitWin,hitWin),2),mean(tsimIndex(missWin,missWin),2),mean(tsimIndex(FAWin,FAWin),2)];
figure,customBoxplot(ensembleStab),box off
set(gca,'TickDir','out','fontsize',14),ylabel('Ensemble Stability')
%% Volume analysis
%%% Fraction of shared L2/3 and L5 neurons within an ensemble
Ensemble = FAEnsemble;
for n = 1:Ensemble.ensembleNum
    temp = Ensemble.rankEnsembles{n};
    l23Member = sum(ismember(temp,l23Idx));
    l5Member = sum(ismember(temp,L5Idx));
    ensembleLayerFrac(n,:) = [l23Member/(l23Member+l5Member),l5Member/(l23Member+l5Member)];
end
figure,pie(mean(ensembleLayerFrac))
%% Ensembles shared across task variables
% Calculate most occurance neurons across ensembles
Ensemble = FAEnsemble;
temp = cellfun(@(x) length(x),Ensemble.rankEnsembles);
tEnsemble = vertcat(Ensemble.rankEnsembles{temp>2});
[uE,~,ix] = unique(tEnsemble);
C = accumarray(ix,1);
[~,idx] = sort(C,'descend');
sharedEFA = [uE(idx),C(idx)]; % Number of reoccuring neurons across ensembles

%%% now we have to calculate reoccurance from within and between conditions
% As an aside we can just loop through the values and find which values are
% members in which conditions then use the frequency of such neurons within
% conditions as weighted values. Ie. If a neuron is shared between many
% ensembles and is consistent across task variables then it is essentially
% TASK INVARIANT meaning it encoding lots of different variables 

sharedEhm = ismember(sharedEhit(:,1),sharedEmiss(:,1)); %shared ensembles between hit and miss
sharedEhf = ismember(sharedEhit(:,1),sharedEFA(:,1)); %shared ensembles between hit and FA
sharedEfm = ismember(sharedEFA(:,1),sharedEmiss(:,1));

sharedEacrossTask = [sum(sharedEhm)/length(sharedEhm),sum(sharedEhf)/length(sharedEhf),sum(sharedEfm)/length(sharedEfm)];
figure,bar(sharedEacrossTask)