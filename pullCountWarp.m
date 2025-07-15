%%
function [warped_lever, warping_functions] = pullCountWarp(pullCount, fs, parameters,Spikes)

%% Downsample pullCount to match calcium data
ca_samples = size(Spikes.PSTH.MIHit.spkRate,2);
ds_pullCount = resample((pullCount*10)', ca_samples, size(pullCount,2))';
ds_pullCount = round(ds_pullCount/10);
[num_trials, num_samples] = size(ds_pullCount);
time_axis = linspace(0,parameters.windowBeforeMI+parameters.windowAfterMI,num_samples)*1000; % ms
% Convert windowBeforeMI to milliseconds
MOVEMENT_OFFSET = parameters.windowBeforeMI * 1000; % Convert seconds to ms

% Create fixed output time grid
output_time = time_axis;
warped_lever = nan(size(ds_pullCount));
warping_functions = cell(num_trials, 1);

%% 1. Detect the ACTUAL pull events (not motion initiation)
event_times = cell(num_trials, 1);
valid_trials = [];

for trial = 1:num_trials
    % Find motion initiation (0 appears)
    motion_start = find(ds_pullCount(trial,:) == 0, 1, 'first');
    if isempty(motion_start), continue; end
    
    % Find ALL pull transitions (increases in pull count)
    pull_diff = diff(ds_pullCount(trial, motion_start:end));
    pull_increases = find(pull_diff > 0) + motion_start;
    
    % We need exactly 3 pulls: 0->1, 1->2, 2->3
    if numel(pull_increases) >= 3
        % Store the three ACTUAL pull events (not motion initiation)
        event_times{trial} = time_axis(pull_increases(1:3)); % [pull1, pull2, pull3]
        valid_trials = [valid_trials, trial];
    end
end

% Define target times for the three ACTUAL pulls
pl = ((mean(Spikes.PSTH.MIHit.pl)-45)/fs)*1000;
%pl = pl-pl(1);
TARGET_PULL1 = pl(1);    % First pull (0->1) at movement offset  
TARGET_PULL2 = pl(2);  % Second pull (1->2) at +500ms
TARGET_PULL3 = pl(3); % Third pull (2->3) at +1000ms

fprintf('Alignment targets with %.1fs movement offset:\n', parameters.windowBeforeMI);
fprintf('Pull 1 (0->1): %.0fms, Pull 2 (1->2): %.0fms, Pull 3 (2->3): %.0fms\n', ...
    MOVEMENT_OFFSET + TARGET_PULL1, MOVEMENT_OFFSET + TARGET_PULL2, MOVEMENT_OFFSET + TARGET_PULL3);

%% 2. Create warping to align the ACTUAL pull events
for i = 1:numel(valid_trials)
    trial = valid_trials(i);
    events = event_times{trial}; % [pull1_time, pull2_time, pull3_time]
    
    % Create inverse warping: for each output time, find corresponding input time
    sample_times = zeros(size(output_time));
    
    for t = 1:length(output_time)
        current_time = output_time(t);
        
        % Define target times with movement offset
        target_pull1 = MOVEMENT_OFFSET + TARGET_PULL1;  % When 0->1 should occur
        target_pull2 = MOVEMENT_OFFSET + TARGET_PULL2;  % When 1->2 should occur  
        target_pull3 = MOVEMENT_OFFSET + TARGET_PULL3;  % When 2->3 should occur
        
        if current_time < target_pull1
            % Before first pull: map to time before original first pull
            time_before = current_time - target_pull1;
            sample_times(t) = events(1) + time_before;
            
        elseif current_time < target_pull2
            % Pull1 to pull2: map [target_pull1, target_pull2] -> [events(1), events(2)]
            progress = (current_time - target_pull1) / (TARGET_PULL2 - TARGET_PULL1);
            sample_times(t) = events(1) + progress * (events(2) - events(1));
            
        elseif current_time < target_pull3
            % Pull2 to pull3: map [target_pull2, target_pull3] -> [events(2), events(3)]
            progress = (current_time - target_pull2) / (TARGET_PULL3 - TARGET_PULL2);
            sample_times(t) = events(2) + progress * (events(3) - events(2));
            
        else
            % After pull3: maintain last segment slope
            excess_time = current_time - target_pull3;
            slope = (events(3) - events(2)) / (TARGET_PULL3 - TARGET_PULL2);
            sample_times(t) = events(3) + excess_time * slope;
        end
    end
    
    % Sample original data at computed time points
    warped_lever(trial,:) = interp1(time_axis, ds_pullCount(trial,:), ...
                                   sample_times, 'nearest', 'extrap');
    
    % Store warping function
    %warping_functions{trial} = @(data) interp1(time_axis, data, sample_times, 'linear', 'extrap');
    % With this:
    %warping_functions{trial} = @(t, data) interp1(t, data, sample_times, 'linear', 'extrap');
    % With this updated version that accepts time parameter:
    warping_functions{trial} = @(input_time, data) flexible_warp(time_axis, sample_times, input_time, data);
end

%% Handle invalid trials
invalid_trials = setdiff(1:num_trials, valid_trials);
for trial = invalid_trials
    warped_lever(trial,:) = ds_pullCount(trial,:);
    warping_functions{trial} = @(x) x;
end

%% Enhanced visualization
figure('Color','white','Position',[100,100,700,500]);

% Original data
subplot(3,1,1);
imagesc(ds_pullCount(valid_trials,:));
title('Original Lever Data');
xlabel('Time (samples)'); ylabel('Trial');
colorbar; colormap('turbo');

% Warped data with alignment lines
subplot(3,1,2);
imagesc(warped_lever(valid_trials,:));
title('Perfectly Aligned Lever Data (ALL 3 Pull Events)');
xlabel('Time (samples)'); ylabel('Trial');
colorbar; colormap('turbo');

% Add vertical lines at target alignment points
hold on;
target_samples = [MOVEMENT_OFFSET + TARGET_PULL1, MOVEMENT_OFFSET + TARGET_PULL2, MOVEMENT_OFFSET + TARGET_PULL3] * fs / 1000;
for i = 1:length(target_samples)
    xline(target_samples(i), 'w--', 'LineWidth', 2);
end

% Individual traces with alignment markers
subplot(3,1,3);
hold on;
for i = 1:min(10, numel(valid_trials))
    trial = valid_trials(i);
    plot(warped_lever(trial,:) + 0.2*i, 'LineWidth', 1.5);
end
title('Individual Warped Traces - Perfect Alignment of ALL 3 Pull Events');
xlabel('Time (samples)'); ylabel('Pull Count (Offset)');

% Add vertical lines at target times
for i = 1:length(target_samples)
    xline(target_samples(i), 'k--', 'LineWidth', 1.5, 'Alpha', 0.7);
end
grid on;

% Add text labels for the alignment targets
text(target_samples(1), max(ylim)*0.9, 'Pull 1 (0→1)', 'Rotation', 90);
text(target_samples(2), max(ylim)*0.9, 'Pull 2 (1→2)', 'Rotation', 90);
text(target_samples(3), max(ylim)*0.9, 'Pull 3 (2→3)', 'Rotation', 90);

fprintf('Warping complete. ALL 3 pull events are now perfectly aligned:\n');
fprintf('- Pull 1 (0→1) at sample %.0f\n', target_samples(1));
fprintf('- Pull 2 (1→2) at sample %.0f\n', target_samples(2)); 
fprintf('- Pull 3 (2→3) at sample %.0f\n', target_samples(3));

end

function warped_data = flexible_warp(lever_time_axis, lever_sample_times, neural_time_axis, neural_data)
% FLEXIBLE_WARP Warps neural data using lever-derived time mapping
% Inputs:
%   lever_time_axis: Original time axis from lever data (ms)
%   lever_sample_times: Warped time points from lever alignment (ms)  
%   neural_time_axis: Time axis for neural data (ms)
%   neural_data: Neural trace to be warped
% Output:
%   warped_data: Warped neural data

    % Create the time mapping function from lever warping
    % This maps: original_time -> warped_time
    time_mapping_func = @(t) interp1(lever_time_axis, lever_sample_times, t, 'linear', 'extrap');
    
    % Apply time mapping to neural data's time axis
    warped_neural_times = time_mapping_func(neural_time_axis);
    
    % Add reflective padding to neural data to handle extrapolation
    padding_length = round(0.1 * length(neural_data));
    
    % Reflective padding
    if length(neural_data) > padding_length
        left_pad = neural_data(padding_length:-1:1);
        right_pad = neural_data(end:-1:end-padding_length+1);
    else
        left_pad = repmat(neural_data(1), 1, padding_length);
        right_pad = repmat(neural_data(end), 1, padding_length);
    end
    
    padded_neural_data = [left_pad, neural_data, right_pad];
    
    % Create padded time axis
    dt = mean(diff(neural_time_axis));
    left_time = neural_time_axis(1) - (length(left_pad):-1:1) * dt;
    right_time = neural_time_axis(end) + (1:length(right_pad)) * dt;
    padded_neural_time = [left_time, neural_time_axis, right_time];
    
    % Constrain warped times to avoid excessive extrapolation
    time_range = max(padded_neural_time) - min(padded_neural_time);
    extrapolation_limit = 0.05 * time_range; % 5% extrapolation limit
    
    warped_neural_times_safe = max(min(warped_neural_times, ...
                                     max(padded_neural_time) - extrapolation_limit), ...
                                 min(padded_neural_time) + extrapolation_limit);
    
    % Apply final interpolation
    warped_data = interp1(padded_neural_time, padded_neural_data, ...
                         warped_neural_times_safe, 'linear', 'extrap');
end