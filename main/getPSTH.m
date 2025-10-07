%% make leverPSTH for self-initiated trials
% Pass through deconvolved calcium data to extract trial aligned data
% across all neurons
function Spikes = getPSTH(spks,Behaviour)
spksTime = linspace(0,size(spks,2)/30.048,size(spks,2));
spks = smoothdata(spks,2,'gaussian',10);
spks = spks-mean(spks,'all');
spks = spks/std(spks,[],'all');
%%% REWARD ALIGNED
neurons = [];
% Here we fix the arduino drift if neccessary
arduinoDrift = (spksTime(end)/(Behaviour.B(end,2)/10^6));
for n = 1:length(Behaviour.hitTrace)
    targetTime = (Behaviour.hitTrace(n).t0-Behaviour.parameters.windowBeforePull)*arduinoDrift;
    % assign index with nearest time index
    st = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    st = max(1, min(length(spksTime), st));

    targetTime = (Behaviour.hitTrace(n).t0+Behaviour.parameters.windowAfterPull)*arduinoDrift;
    % assign index with nearest time index
    sp = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    sp = max(1, min(length(spksTime), sp));

    neurons{n} = spks(:,st:sp); % grab all neurons for the trial
    for sQN = 1:(length(unique(Behaviour.hitTrace(n).pullCount))-1) % substract 1 because we dont want zeros
        [~,idx] = min(abs(Behaviour.hitTrace(n).rawtime-Behaviour.hitTrace(n).t0));
        tT = find(Behaviour.hitTrace(n).pullCount==sQN);
        tT(tT<idx) = [];
        tT = tT(1);
        % Index of pull
        targetTime = (Behaviour.hitTrace(n).rawtime(tT))*arduinoDrift;
        pl = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
        % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
        pl = max(1, min(length(spksTime), pl));
        pl = pl-st;
        Spikes.PSTH.hit.pl(n,sQN) = pl;
    end
end

trialLen = mode(cellfun(@(x) size(x,2),neurons));
neurons = cellfun(@(x) x(:,1:trialLen),neurons,'UniformOutput',false);

% Get the number of neurons from the first trial
numNeurons = size(neurons{1}, 1);

% Get the time points per trial (assuming all trials have the same length)
timePerTrial = size(neurons{1}, 2);

% Concatenate all trial matrices along the second dimension (time)
% The result will be numNeurons x (sum of time points across trials)
reshaped_output = cell2mat(neurons);
trials = [];
for n = 1:numNeurons
    trials{n} = reshape(reshaped_output(n,:),timePerTrial,[])';
end

output = make_nice_mean_raster(trials,1,0);
Spikes.PSTH.hit.spks = neurons;
Spikes.PSTH.hit.spkRate = output;
Spikes.PSTH.hit.trialRate = make_nice_mean_raster(neurons,1,0);

neurons = [];
for n = 1:length(Behaviour.missTrace)
    targetTime = (Behaviour.missTrace(n).t0-Behaviour.parameters.windowBeforePull)*arduinoDrift;
    % assign index with nearest time index
    st = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    st = max(1, min(length(spksTime), st));

    targetTime = (Behaviour.missTrace(n).t0+Behaviour.parameters.windowAfterPull)*arduinoDrift;
    % assign index with nearest time index
    sp = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    sp = max(1, min(length(spksTime), sp));
    neurons{n} = spks(:,st:sp);
end

trialLen = mode(cellfun(@(x) size(x,2),neurons));
neurons = cellfun(@(x) x(:,1:trialLen),neurons,'UniformOutput',false);
% Get the number of neurons from the first trial
numNeurons = size(neurons{1}, 1);

% Get the time points per trial (assuming all trials have the same length)
timePerTrial = size(neurons{1}, 2);

% Concatenate all trial matrices along the second dimension (time)
% The result will be numNeurons x (sum of time points across trials)
reshaped_output = cell2mat(neurons);
trials = [];
for n = 1:numNeurons
    trials{n} = reshape(reshaped_output(n,:),timePerTrial,[])';
end

output = make_nice_mean_raster(trials,1,0);
Spikes.PSTH.miss.spks = neurons;
Spikes.PSTH.miss.spkRate = output;
Spikes.PSTH.miss.trialRate = make_nice_mean_raster(neurons,1,0);

%%% MOVEMENT ALIGNED
neurons = [];
for n = 1:length(Behaviour.MIHitTrace)
    targetTime = (Behaviour.MIHitTrace(n).t0-Behaviour.parameters.windowBeforeMI)*arduinoDrift;
    % assign index with nearest time index
    st = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    st = max(1, min(length(spksTime), st));

    targetTime = (Behaviour.MIHitTrace(n).t0+Behaviour.parameters.windowAfterMI)*arduinoDrift;
    % assign index with nearest time index
    sp = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
    % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
    sp = max(1, min(length(spksTime), sp));
    neurons{n} = spks(:,st:sp);
    for sQN = 1:(length(unique(Behaviour.hitTrace(n).cleanedPullCounts))-1) % substract 1 because we dont want zeros
        [~,idx] = min(abs(Behaviour.MIHitTrace(n).rawtime-Behaviour.MIHitTrace(n).t0));
        tT = find(Behaviour.MIHitTrace(n).cleanedPullCounts==sQN);
        tT(tT<idx) = [];
        tT = tT(1);
        % Index of pull
        targetTime = (Behaviour.MIHitTrace(n).rawtime(tT))*arduinoDrift;
        pl = round(interp1(spksTime, 1:length(spksTime), targetTime, 'nearest', 'extrap'));
        % Ensure indices are within bounds (interp1 with 'extrap' might produce out-of-bounds indices)
        pl = max(1, min(length(spksTime), pl));
        pl = pl-st;
        Spikes.PSTH.MIHit.pl(n,sQN) = pl;
    end
end

trialLen = mode(cellfun(@(x) size(x,2),neurons));
neurons = cellfun(@(x) x(:,1:trialLen),neurons,'UniformOutput',false);
% Get the number of neurons from the first trial
numNeurons = size(neurons{1}, 1);

% Get the time points per trial (assuming all trials have the same length)
timePerTrial = size(neurons{1}, 2);

% Concatenate all trial matrices along the second dimension (time)
% The result will be numNeurons x (sum of time points across trials)
reshaped_output = cell2mat(neurons);
trials = [];
for n = 1:numNeurons
    trials{n} = reshape(reshaped_output(n,:),timePerTrial,[])';
end

output = make_nice_mean_raster(trials,20,0);
Spikes.PSTH.MIHit.spks = neurons;
Spikes.PSTH.MIHit.spkRate = output;
Spikes.PSTH.MIHit.trialRate = make_nice_mean_raster(neurons,1,0);
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


