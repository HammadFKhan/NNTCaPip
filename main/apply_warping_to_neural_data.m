function warped_neural_cell = apply_warping_to_neural_data(neural_cell, warping_functions, neural_fs)
% Inputs:
%   neural_cell: {trial} cell array, each [neurons x time]
%   warping_functions: {trial} cell array of function handles
%   neural_fs: Sampling rate of neural data (Hz)

num_trials = length(neural_cell);
warped_neural_cell = cell(size(neural_cell));

for trial = 1:num_trials
    if isempty(warping_functions{trial})
        warped_neural_cell{trial} = neural_cell{trial};
        continue;
    end
    
    trial_data = neural_cell{trial};
    [num_neurons, num_timepoints] = size(trial_data);
    
    % Create neural-specific time axis (milliseconds)
    neural_time_axis = (0:num_timepoints-1) / neural_fs * 1000; 
    
    warped_trial_data = [];
    
    for neuron = 1:num_neurons
        neural_trace = trial_data(neuron, :);
        
        % Pass BOTH time vector and data to warping function
        try
        warped_trace = warping_functions{trial}(neural_time_axis, neural_trace);
        warped_trial_data = [warped_trial_data;warped_trace];
        catch ME
            disp('Skipped trial')
            break
        end
    end
    
    warped_neural_cell{trial} = warped_trial_data;
end
