clear all
addpath(genpath('main'));
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
%     
    % Fix centroids
    ROIcentroid = [];
    for i = 1:length(ROI)
        if isempty(ROI{i}), ROI{i} = {[0 0]};end
        blah = vertcat(ROI{i}{:});
        ROIcentroid(i,:) = floor(mean(blah,1));
    end
    set(0,'DefaultFigureWindowStyle','normal')
    addpath(genpath('main'));
    addpath(genpath('Pipelines'));
    std_threshold = 7;
    static_threshold = .01;
    Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
    %Excude inactive cells
    % numSpikes = sum(Spikes,2);
    % keepSpikes = find(numSpikes>(.01*mean(numSpikes)));
    % Spikes = Spikes(keepSpikes,:);
    [coactive_cells,detected_spikes] = coactive_index(Spikes,5000);
    cell_count = length(ROI);
    time = time_adjust(size(DeltaFoverF,2),30.048);
    for i = 1:size(DeltaFoverF,1)
        try
        calcium_avg{i} = STA(DeltaFoverF(i,:),2,120);%std, window (frames)
        catch ME
            continue;
        end
    end
    
    % Perform shuffling and pairwise if data is small enough
    if size(DeltaFoverF,2)<2000
        %     Spikes_shuffled = tempShuffle(Spikes,1000);
        %     Event_shuffled = spatialShuffle(Spikes,1000);
        %     surrogate = 10;
        %     Total_shuffled = allShuffle(Spikes,1000);
        %     [shufcoactive_cells,detected_spikes] = coactive_index(Spikes_shuffled,length(Spikes_shuffled));
        %     shuff_corr = correlation_dice(Event_shuffled);
        %     [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,bin);
        %     shufsim_index = shufsim_index-mean(mean(shufsim_index,2));
        %     factorCorrection = 100*floor(size(Spikes,2)/100); % Correct for frame size aquisition
        %     [vectorized,sim_index] = cosine_similarity(Spikes(:,1:factorCorrection),50);
        corr = correlation_dice(Spikes);
        Connected_ROI = Connectivity_dice(corr,0.1);
        [NumActiveNodes,NodeList,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,...
            ActivityCentroid,ActivityCentroidVariance]...
            = Network_Analysis(ROIcentroid,Connected_ROI);
    end
    % Pairwise Velocity Analysis
    % velocityPairwise(VR_data,Spikes)
    % Ensemble Analysis
    % figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
    factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
    Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);
    
    % Ensemble stats
    Ensemble = ensembleMetric(Ensemble,AverageImage,ROIcentroid);
    Ensemble = ensembleStat(Ensemble);
    close all
    
    if ~exist([file(fileNum).folder '\Dendritic'],'dir')
        mkdir([file(fileNum).folder '\Dendritic']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\output'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\Dendritic'],[file_name,'.mat']);
    save(savepath,'files', 'Ensemble','Spikes','ROI', 'ROIcentroid' ,'DeltaFoverF');
    clearvars -except file numFile fileNum filetype foldername
end