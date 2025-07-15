function warped_lever = warpLever(lever_data, warping_functions, fs)
% Inputs:
%   neural_cell: {trial} cell array, each [neurons x time]
%   warping_functions: {trial} cell array of function handles
%   neural_fs: Sampling rate of neural data (Hz)


[num_trials, num_timepoints] = size(lever_data);

% Create neural-specific time axis (milliseconds)
lever_time_axis = (0:num_timepoints-1) / fs * 1000;
warped_lever = nan(size(lever_data));
for trial = 1:num_trials
    lever_trace = lever_data(trial, :);
    % Pass BOTH time vector and data to warping function
    try
        warped_trace = warping_functions{trial}(lever_time_axis, lever_trace);
        warped_lever(trial,:) = warped_trace;
    catch ME
        disp('Skipped trial')
        continue
    end
end
