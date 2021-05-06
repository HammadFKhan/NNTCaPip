%% Data Analysis for coactive cells across Trails
pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
directory = dir(pathname);
L = length(directory);

for i = 3:L
     fullfile = strcat(directory(i).folder,'\',directory(i).name);
     load(fullfile);
     k = i-2;
    std_threshold = 2.5;
    static_threshold = 0.2;
    Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
    Spikes_shuffled = tempShuffle(Spikes,num_images,cell_count);
    [coactive_cells,detected_spikes] = coactive_index(Spikes_shuffled);

    list(k,1) = convertCharsToStrings(directory(i).name);
    list1(k,1) = mean(coactive_cells);
    list1(k,2) = max(coactive_cells);
%     out.coactive_cells_run = movmean(coactive_cells,100);coactive_cells_run_avg = mean(coactive_cells_run);
end
%% Mean Velocity across trials
pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
directory = dir(pathname);
L = length(directory);
for ii = 3:L
    k = ii-2;
    fullfile = strcat(directory(ii).folder,'\',directory(ii).name);
    encoder_data = convert_encoder(fullfile);
    nameList(k,1) = convertCharsToStrings(directory(ii).name);
    vList(k,1) = mean(encoder_data.ang_velocity);
    if isnan(vList(k,1))
        vList(k,1) = 0;
    end
end
vMean = mean(vList);