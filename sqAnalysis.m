%% Pipeline to analyze lever pull without cue response
parameters.experiment = 'self'; % self - internally generated, cue - cue initiated
parameters.opto = 0; % 1 - opto ON , 0 - opto OFF
parameters.cool = 0; % No Cool 
parameters.windowBeforePull = 3; % in seconds
parameters.windowAfterPull = 2; % in seconds
parameters.windowBeforeCue = 1.5; % in seconds
parameters.windowAfterCue = 1.5; % in seconds
parameters.windowBeforeMI = 1.5; % in seconds 
parameters.windowAfterMI = 3.5; % in seconds 
parameters.delay = 0.5; %reward delay
parameters.Fs = 1000; % Eventual downsampled data
parameters.ts = 1/parameters.Fs;
parameters.rows = 64;
parameters.cols = 1;

%bFile = 'D:\SQLever\Seq\Mouse2Day19\Rbp4Mouse2Day192025_05_28_12.02.PM.csv';
bFile = 'D:\SQLever\Seq\Mouse2Day28Sq\SqMouse2Day282025_06_04_01.57.PM.csv';
addpath(genpath('C:\Users\khan332\Documents\GitHub\NNTEphysPip\Main\BehaviourAnalysis'))
%[Behaviour] = readLever(parameters,[],bFile);
[Behaviour] = readLeverSq(parameters,[],bFile);
Behaviour.parameters = parameters;
%% Behaviour
% Reward Aligned
figure;
time = linspace(-parameters.windowBeforePull, parameters.windowAfterPull,length(Behaviour.hitTrace(1).rawtrace));
dat = arrayfun(@(x) mean(x.rawtrace,2),Behaviour.hitTrace,'UniformOutput',false);
leverTraces = horzcat(dat{:})';

plot(time,smoothdata(leverTraces,2,'movmean',10),'color',[0.5 0.5 0.5 0.25]),hold on
xline(-0.5,'r--','Sq Complete')
xline(0.0,'r--','Reward')

plot(time,smoothdata(mean(leverTraces,1),'movmean',10),'r')

allPulls = arrayfun(@(x) x.pullCount, Behaviour.hitTrace, 'UniformOutput', false);

% Determine the correct size (number of rows) from the first array
correctNumRows = size(allPulls{1}, 1);

% Find which arrays have the correct number of rows
validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

% Keep only valid arrays
validPulls = allPulls(validIdx);

% Horizontally concatenate and transpose as you did
allPulls = horzcat(validPulls{:})';
pullCounts = allPulls;
% [cleanedpullCounts, hasTimeout] = cleanTimeoutSequences(pullCounts,Behaviour);
figure,
imagesc([-parameters.windowBeforePull*1000 parameters.windowAfterPull*1000], [1 size(pullCounts,1)], pullCounts);
% Assume pullCounts is [trials x time] and timeAxis is the time vector (e.g., -1500:bin:1500)
[numTrials, numBins] = size(pullCounts);
timeAxis = linspace(-parameters.windowBeforeMI , parameters.windowAfterMI, numBins); % adjust as needed

durations = nan(numTrials,1);

for t = 1:numTrials
    % Only consider pulls before reward (time < 0)
    pulls = find(cleanedpullCounts(t,:) > 0 & timeAxis > 0);
    if ~isempty(pulls)
        durations(t) = timeAxis(pulls(end)) - timeAxis(pulls(1));
    end
end

% Sort by duration (ascending: fastest to slowest)
[~, sortIdx] = sort(durations, 'ascend', 'MissingPlacement','last'); % NaNs (no pulls) go to bottom

% Reorder pullCounts for plotting
sortedPullCounts = cleanedpullCounts(sortIdx, :);

% Plot
figure;
imagesc(timeAxis, 1:numTrials, sortedPullCounts);
colormap('parula'); % or your preferred colormap
colorbar;
xlabel('Time from reward (ms)');
ylabel('Trials (fastest to slowest)');
title('Lever Pulls Sorted by Sequence Duration');
%% first movement aligned
% 
figure;
time = linspace(-parameters.windowBeforeMI, parameters.windowAfterMI,length(Behaviour.MIHitTrace(1).rawtrace));
leverTraces = horzcat(Behaviour.MIHitTrace.rawtrace)';

plot(time,smoothdata(leverTraces,2,'movmean',10),'color',[0.5 0.5 0.5 0.25]),hold on
xline(0.0,'r--','Movement')
xline(mean(vertcat(Behaviour.MIHitTrace.rewardtime)),'r--','Reward')

plot(time,smoothdata(mean(leverTraces,1),'movmean',10),'r')

allPulls = arrayfun(@(x) x.pullCount, Behaviour.MIHitTrace, 'UniformOutput', false);

% Determine the correct size (number of rows) from the first array
correctNumRows = size(allPulls{1}, 1);

% Find which arrays have the correct number of rows
validIdx = cellfun(@(c) size(c,1) == correctNumRows, allPulls);

% Keep only valid arrays
validPulls = allPulls(validIdx);

% Horizontally concatenate and transpose as you did
allPulls = horzcat(validPulls{:})';
pullCounts = allPulls;
[cleanedpullCounts, Behaviour] = cleanTimeoutSequences(pullCounts,Behaviour);
figure,
imagesc([-parameters.windowBeforeMI parameters.windowAfterMI], [1 size(pullCounts,1)], (cleanedpullCounts));

%% Load in calcium data and clean up
%load('D:\SQLever\Seq\Mouse2Day19\denoised_output\suite2p\plane0\Fall.mat')
%load('D:\SQLever\Seq\Mouse2Day19\denoised_output\DataFolderIs_chunks_202506021548_ModelFolderIs_soma_best_model\E_10_Iter_6048\suite2p\plane0\Fall.mat')
load('D:\SQLever\Seq\Mouse2Day28Sq\chunks\DataFolderIs_chunks_202506081421_ModelFolderIs_soma_best_model\E_20_Iter_6064\suite2p\plane0\Fall.mat')
ROIcentroid = cellfun(@(x) x.med, stat,'UniformOutput',false);
ROIcentroid = vertcat(ROIcentroid{:});

%%% Clean up soma
F1 = F(iscell(:,1)==1,:);
spks1 = spks(iscell(:,1)==1,:);
skew = cellfun(@(x) x.skew, stat);
skew = skew(iscell(:,1)==1);
keepId = (skew>1.5)';
Fnew = F1(keepId,:);
Fneun = Fneu(keepId,:);
spks1 = spks1(keepId,:);
ROIcentroid = ROIcentroid(keepId,:);
%%% Run convolved traces after loading in suite2P data
addpath(genpath('main'))
addpath(genpath('CaImAn-MATLAB'))
[dF,F_detrended] = convolveFluroescence(Fnew,Fneun,0); % no convolve just detrend
% Segment neurons into individual trial data
Spikes = getPSTH(spks1,Behaviour);
% Spikes = getSqPSTH(spks1,Behaviour,Spikes);
Spikes = getSqPSTH2(spks1,Behaviour,Spikes);


%% Curate data for warping analysis in Python
SqSpikes = zeros(length(Spikes.PSTH.hit.spks),size(Spikes.PSTH.hit.spks{1},2),size(Spikes.PSTH.hit.spks{1},1));
for n = 1:size(SqSpikes,1)
    SqSpikes(n,:,:) = Spikes.PSTH.hit.spks{n}';
end
pullIndex = Spikes.PSTH.MIHit.pl;
rmIndex = find(sum(pullIndex,2)==0);
pullIndex(rmIndex,:) = [];
SqSpikes(rmIndex,:,:) = [];
[fpath,name,exts] = fileparts(bFile);
sessionName = [fpath,'\','SqCa.mat'];
save(sessionName,"SqSpikes","pullIndex","fpath");
disp('Data Saved for warping')
%% Plot out PSTH of neurons
time = linspace(-parameters.windowBeforePull, parameters.windowAfterPull,size(Spikes.PSTH.hit.spkRate,2));
figure,plot(time,Spikes.PSTH.hit.trialRate,'color',[0.0 0.0 0.5 0.25]),hold on
xline(0.0,'k','Sq end')
figure,plot(time,Spikes.PSTH.miss.trialRate,'color',[0.5 0.0 0.0 0.25]),hold on
time = linspace(-parameters.windowBeforeMI, parameters.windowAfterMI,size(Spikes.PSTH.hit.spkRate,2));
figure,plot(time,Spikes.PSTH.MIHit.trialRate,'color',[0.0 0.5 0.0 0.25]),hold on
xline(0.0,'r--','Movement')
xline(mean(vertcat(Behaviour.MIHitTrace.rewardtime)),'r--','Reward')
plot(time,mean(Spikes.PSTH.MIHit.trialRate),'k')
%figure,plot(time,mean(Spikes.PSTH.hit.spkRate,1),'k')
%%
time = linspace(-parameters.windowBeforePull, parameters.windowAfterPull,size(Spikes.PSTH.hit.spkRate,2));
figure,plot(time,Spikes.PSTH.hit.spkRate,'color',[0.0 0.0 0.5 0.25]),hold on
figure,plot(time,Spikes.PSTH.miss.spkRate,'color',[0.5 0.0 0.0 0.25]),hold on
time = linspace(-parameters.windowBeforeMI, parameters.windowAfterMI,size(Spikes.PSTH.hit.spkRate,2));
figure,plot(time,Spikes.PSTH.MIHit.spkRate,'color',[0.0 0.5 0.0 0.25]),hold on
%% 
id = 2;
figure,
exampleSpk = cellfun(@(x) x(id,:), Spikes.PSTH.hit.spks,'UniformOutput',false);
exampleSpk = vertcat(exampleSpk{:});
subplot(1,3,1),imagesc(time,1:size(exampleSpk,1),exampleSpk),caxis([0 1]),axis square
exampleSpk = cellfun(@(x) x(id,:), Spikes.PSTH.miss.spks,'UniformOutput',false);
exampleSpk = vertcat(exampleSpk{:});
subplot(1,3,2),imagesc(time,1:size(exampleSpk,1),exampleSpk),caxis([0 1]),axis square
exampleSpk = cellfun(@(x) x(id,:), Spikes.PSTH.MIHit.spks,'UniformOutput',false);
exampleSpk = vertcat(exampleSpk{:});
subplot(1,3,3),imagesc(time,1:size(exampleSpk,1),exampleSpk),caxis([0 1]),axis square

%%
cmap = cmocean('balance');
time = linspace(-Spikes.SQ_PSTH(1).win,Spikes.SQ_PSTH(1).win,size(Spikes.SQ_PSTH(1).MIHit.trialRate,2));
figure,
for n = 1:3
    subplot(1,3,n),imagesc(time,1:334,zscore(Spikes.SQ_PSTH(n).MIHit.spkRate')'),colormap(cmap),caxis([-2.54 2.54]),axis square
end

figure,
for n = 1:3
    subplot(1,3,n),imagesc(time,1:91,(Spikes.SQ_PSTH(n).MIHit.trialRate)),colormap(jet),axis square
end


figure,hold on
for n = 1:3
    subplot(1,3,n),plot(time,mean(Spikes.SQ_PSTH(n).MIHit.trialRate)')
    xline(0,'k','movement')
    if n == 3
        xline(15,'r','reward')
    end
end
legend('Pull 1','Pull 2','Pull 3')



figure,hold on
for n = 1:3
    subplot(1,3,n),plot(time,(Spikes.SQ_PSTH(n).MIHit.spkRate)')
end
legend('Pull 1','Pull 2','Pull 3')
%%
% Assuming:
% neural_cell: {trial1, trial2, ...} where trial1 = [neurons x time] matrix
% warping_functions: from perfect_pull_alignment_fixed()
neural_cell = Spikes.PSTH.hit.spks;
% Apply warping
warped_neural_cell = apply_warping_to_neural_data(neural_cell, warping_functions, 30);
emptyCells = cellfun(@isempty, warped_neural_cell);
warped_neural_cell(emptyCells) = [];
%% Verify with a sample neuron and trial
dat = cellfun(@(x) x(201,:),neural_cell,'UniformOutput',false);
dat = vertcat(dat{:});

figure,plot(dat')
clim([0 500])
%% Load in warped data from affinewarp
load('D:\SQLever\Seq\Mouse2Day28Sq\warpedSpks.mat')
disp('Loaded warp spikes data')

figure,
for n = 1:15
    subplot(5,3,n),imagesc(squeeze(warpedSpks.leverWarpedSpks_sorted(:,:,n))),colorbar
end
warpedSpks = getSqPSTHwarped(warpedSpks,Behaviour,Spikes);
%%
cmap = cmocean('balance');
time = linspace(-warpedSpks.SQ_PSTH(1).win,warpedSpks.SQ_PSTH(1).win,size(warpedSpks.SQ_PSTH(1).MIHit.trialRate,2));
figure,
for n = 1:3
    subplot(1,3,n),imagesc(time,1:334,zscore(warpedSpks.SQ_PSTH(n).MIHit.spkRate')'),colormap(cmap),caxis([-2.54 2.54]),axis square
end

figure,
imagesc((squeeze(warpedSpks.leverWarpedSpks_sorted(:,:,2)))),colormap(cmap),axis square


figure
for n = 1:3
    plot(time,mean(warpedSpks.SQ_PSTH(n).MIHit.trialRate)'),hold on,axis square
end
xline(0,'k','movement')
xline(15,'r','reward')



figure,hold on
for n = 1:3
    subplot(3,1,n),plot(time,(warpedSpks.SQ_PSTH(n).MIHit.spkRate)'),hold on
end
legend('Pull 1','Pull 2','Pull 3')
%% Plot for full trial
time = linspace(-parameters.windowBeforeMI, parameters.windowAfterMI,size(Spikes.PSTH.MIHit.spkRate,2));
figure,subplot(3,1,1),plot(time,Spikes.PSTH.MIHit.trialRate,'color',[0.0 0.0 0.5 0.25]),hold on
xline(0.0,'k','Sq Start')
subplot(3,1,2),plot(time,squeeze(mean(warpedSpks.leverWarpedSpks,3)),'color',[0.0 0.0 0.5 0.25]),hold on
xline(0.0,'k','Sq Start')
subplot(3,1,3),plot(time,squeeze(mean(warpedSpks.warpedSpks,3)),'color',[0.0 0.0 0.5 0.25]),hold on
xline(0.0,'k','Sq Start')
%%
pl = (mean(Spikes.PSTH.MIHit.pl)-45)/30;
pl = pl;
time = linspace(-parameters.windowBeforeMI, parameters.windowAfterMI,size(Spikes.PSTH.MIHit.spkRate,2));
figure,subplot(3,1,1),plot(time,mean(Spikes.PSTH.MIHit.trialRate),'color',[0.0 0.0 0.5 1]),hold on
xline(0.0,'k','Sq Start')
for n = 1:3
    xline(pl(n),'k')
end
subplot(3,1,2),plot(time,squeeze(mean(warpedSpks.leverWarpedSpks,[1,3])),'color',[0.0 0.0 0.5 1]),hold on
xline(0.0,'k','Sq Start')
for n = 1:3
    xline(pl(n),'k')
end
subplot(3,1,3),plot(time,squeeze(mean(warpedSpks.warpedSpks,[1,3])),'color',[0.0 0.0 0.5 1]),hold on
xline(0.0,'k','Sq Start')
for n = 1:3
    xline(pl(n),'k')
end

%% Behavior and spiking activity response

% Example usage:
% [warped_pullCount, warping_functions] = warp_lever_sequences(pullCount, 1000);
% visualize_warping_results(pullCount, warped_pullCount, valid_trials);
% 
% % Apply same warping to neural data
% warped_neural_data = nan(size(neural_data));
% for trial = 1:size(neural_data,1)
%     if ~isempty(warping_functions{trial})
%         warped_neural_data(trial,:) = warping_functions{trial}(neural_data(trial,:));
%     end
% end
[warped_pullCount, warping_functions] = pullCountWarp(cleanedpullCounts, 30, parameters,Spikes);
warped_lever = warpLever(leverTraces, warping_functions, 100);
%% Beautified Warped Lever Visualization
% Create time axis in seconds
[num_trials, num_samples] = size(warped_lever);
time_axis = (0:num_samples-1) / 100; % 100 Hz sampling rate

% Calculate target pull times
MOVEMENT_OFFSET = parameters.windowBeforeMI; % seconds
TARGET_PULL1 = MOVEMENT_OFFSET + 0;      % First pull
TARGET_PULL2 = MOVEMENT_OFFSET + 0.5;    % Second pull  
TARGET_PULL3 = MOVEMENT_OFFSET + 1.0;    % Third pull

figure('Position', [100, 100, 1400, 900], 'Color', 'white');

% Panel 1: Heatmap view
subplot(2,2,1);
imagesc(time_axis, 1:min(50, num_trials), warped_lever(1:min(50, num_trials), :));
colormap(gca, 'turbo');
colorbar;
title('Warped Lever Traces - Heatmap View', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Trial');

% Add alignment lines
hold on;
xline(TARGET_PULL1, 'w--', 'LineWidth', 2, 'Alpha', 0.8);
xline(TARGET_PULL2, 'w--', 'LineWidth', 2, 'Alpha', 0.8);
xline(TARGET_PULL3, 'w--', 'LineWidth', 2, 'Alpha', 0.8);

% Panel 2: Individual traces (offset)
subplot(2,2,2);
hold on;
colors = turbo(min(20, num_trials));
y_offset = 0;

for trial = 1:min(20, num_trials)
    plot(time_axis, warped_lever(trial, :) + y_offset, ...
         'Color', colors(trial, :), 'LineWidth', 1.2);
    y_offset = y_offset + 20; % Offset each trace
end

% Add alignment lines
xline(TARGET_PULL1, 'k--', 'LineWidth', 2, 'Alpha', 0.7);
xline(TARGET_PULL2, 'k--', 'LineWidth', 2, 'Alpha', 0.7);
xline(TARGET_PULL3, 'k--', 'LineWidth', 2, 'Alpha', 0.7);

title('Individual Warped Traces (First 20 Trials)', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Lever Position (Offset)');

% Panel 3: Average trace with confidence intervals
subplot(2,2,3);
focus_start = TARGET_PULL1 - 0.2;
focus_end = TARGET_PULL3 + 0.3;
focus_idx = (time_axis >= focus_start) & (time_axis <= focus_end);
mean_trace = mean(warped_lever, 1, 'omitnan');
std_trace = std(warped_lever, 0, 1, 'omitnan');
sem_trace = std_trace / sqrt(sum(~isnan(warped_lever(:,1))));

% Plot confidence interval
fill([time_axis(focus_idx), fliplr(time_axis(focus_idx))], ...
     [mean_trace(focus_idx) + sem_trace(focus_idx), fliplr(mean_trace(focus_idx) - sem_trace(focus_idx))], ...
     [0.7 0.7 0.9], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
hold on;

% Plot mean trace
plot(time_axis(focus_idx), mean_trace(focus_idx), 'b-', 'LineWidth', 3);

% Add alignment lines with labels
xline(TARGET_PULL1, 'r--', 'Pull 1', 'LineWidth', 2, 'LabelOrientation', 'horizontal');
xline(TARGET_PULL2, 'r--', 'Pull 2', 'LineWidth', 2, 'LabelOrientation', 'horizontal');  
xline(TARGET_PULL3, 'r--', 'Pull 3', 'LineWidth', 2, 'LabelOrientation', 'horizontal');

title('Average Warped Lever Trace ± SEM', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Lever Position');
legend('SEM', 'Mean', 'Location', 'best');

% Panel 4: Pull event magnification
subplot(2,2,4);
% Focus on the pull sequence period
focus_start = TARGET_PULL1 - 0.2;
focus_end = TARGET_PULL3 + 0.3;
focus_idx = (time_axis >= focus_start) & (time_axis <= focus_end);

hold on;
for trial = 1:min(15, num_trials)
    plot(time_axis(focus_idx), warped_lever(trial, focus_idx), ...
         'Color', [0.6 0.6 0.6 0.4], 'LineWidth', 0.8);
end

% Plot average on top
plot(time_axis(focus_idx), mean_trace(focus_idx), 'r-', 'LineWidth', 3);

% Add alignment lines
xline(TARGET_PULL1, 'k--', 'LineWidth', 2);
xline(TARGET_PULL2, 'k--', 'LineWidth', 2);
xline(TARGET_PULL3, 'k--', 'LineWidth', 2);

% Add text labels
text(TARGET_PULL1, max(ylim)*0.9, 'Pull 1', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(TARGET_PULL2, max(ylim)*0.9, 'Pull 2', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(TARGET_PULL3, max(ylim)*0.9, 'Pull 3', 'HorizontalAlignment', 'center', 'FontWeight', 'bold');

title('Pull Sequence Detail', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Time (s)'); ylabel('Lever Position');

sgtitle('Warped Lever Traces - Perfect Temporal Alignment', 'FontSize', 16, 'FontWeight', 'bold');
%%
figure,
for n = 1:25
subplot(5,5,n),imagesc(warped_neural_cell{n}),caxis([0 500]),axis square
end
%%
figure,
imagesc(warped_neural_cell{10}),caxis([0 50]),axis square

figure,
for n = 1:20
    plot(smoothdata(warped_neural_cell{n},2,'gaussian',10)','color',[0.5 0.5 0.5 0.5]), hold on
end
%% Basic functions
function output = make_nice_mean_raster(spmat,smooth_window,showplot)
%*********** spmat1 and spmat2 are spike matrices of two conditions you wish to compare
%*********** smooth_window ... gaussian smoothing in millisecs
numconds = size(spmat,2);
if (numconds==2)
    colo = [[1,0,0];[0,0,1]];
else
    colo = jet(numconds);
end
for k = 1:numconds
    spud = spmat{k};
    numtrials = size(spud,1);
    smorate = gauss_smooth(sum( spud(1:numtrials,:))/....
        numtrials,smooth_window)*1000;
    if showplot
        plot(smorate,'k'); hold on;
        %                 set(H,'Color',colo(k,:));
    end
    output(k,:) = smorate;

end
end

%**************************************************************
function output = gauss_smooth(input, window)
% Smoothing function:
% output = smooth(input, window)
% "Window" is the total kernel width.
% Input array must be one-dimensional.

input_dims = ndims(input);
input_size = size(input);
if input_dims > 2 | min(input_size) > 1,
    disp('Input array is too large.');
    return
end

if input_size(2) > input_size(1),
    input = input';
    toggle_dims = 1;
else
    toggle_dims = 0;
end

if window/2 ~= round(window/2),
    window = window + 1;
end
halfwin = window/2;

input_length = length(input);
%********* gauss window +/- 1 sigma
x = -halfwin:1:halfwin;
kernel = exp(-x.^2/(window/2)^2);
kernel = kernel/sum(kernel);

padded(halfwin+1:input_length+halfwin) = input;
padded(1:halfwin) = ones(halfwin, 1)*input(1);
padded(length(padded)+1:length(padded)+halfwin) = ones(halfwin, 1)*input(input_length);

output = conv(padded, kernel);
output = output(window:input_length+window-1);

if toggle_dims == 1,
    output = output';
end
end
