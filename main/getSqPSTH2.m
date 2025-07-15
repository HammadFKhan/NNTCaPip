function Spikes = getSqPSTH2(spks,Behaviour,Spikes)
%% make leverPSTH for self-initiated trials
% Pass through deconvolved calcium data to extract trial aligned data
% across all neurons
spksTime = linspace(0,size(spks,2)/30.048,size(spks,2));
% Here we fix the arduino drift if neccessary
%%% MOVEMENT ALIGNED
win = floor(Behaviour.parameters.windowBeforeMI*30.048);
neurons = cell(length(Behaviour.MIHitTrace),(length(unique(Behaviour.MIHitTrace(1).cleanedPullCounts))-1));
% Lets align to average lever pull times
dat = Spikes.PSTH.MIHit.pl;
dat(dat==0) = NaN;
pullIdx = floor(nanmean(dat,1));
for n = 1:length(Spikes.PSTH.MIHit.spks)
    for sQN = 1:size(Spikes.PSTH.MIHit.pl,2) % substract 1 because we dont want zeros
         neurons{n,sQN} = Spikes.PSTH.MIHit.spks{n}(:,(pullIdx(sQN)-win):(pullIdx(sQN)+win));
    end
end
emptyCells = cellfun(@isempty, neurons);
neurons(emptyCells(:,1),:) = [];
trialLen = mode(mode(cellfun(@(x) size(x,2),neurons)));
neurons = cellfun(@(x) x(:,1:trialLen),neurons,'UniformOutput',false);
% Get the number of neurons from the first trial
numNeurons = size(neurons{1}, 1);

% Get the time points per trial (assuming all trials have the same length)
timePerTrial = size(neurons{1}, 2);

% Concatenate all trial matrices along the second dimension (time)
% The result will be numNeurons x (sum of time points across trials)
trials = [];
for SqN = 1:size(neurons,2)
    reshaped_output = cell2mat(neurons(:,SqN)');
    for n = 1:numNeurons
        trials{n,SqN} = reshape(reshaped_output(n,:),timePerTrial,[])';
    end
    output = make_nice_mean_raster(trials(:,SqN)',5,0);
    Spikes.SQ_PSTH(SqN).win = win;
    Spikes.SQ_PSTH(SqN).MIHit.spks = neurons;
    Spikes.SQ_PSTH(SqN).MIHit.spkRate = output;
    Spikes.SQ_PSTH(SqN).MIHit.trialRate = make_nice_mean_raster(neurons(:,SqN)',5,0);
    Spikes.SQ_PSTH(SqN).Sequence = SqN;
end

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


