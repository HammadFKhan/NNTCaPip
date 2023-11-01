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
    % Check for bad components
    if exist('badComponents','var') && ~exist('badComFlag','var')
        [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
            removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
        badComFlag = 1;
    end
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
    %     std_threshold = 6;
    %     static_threshold = .00;
    %     Spikes = Spike_Detector_Single((dDeltaFoverF),std_threshold,static_threshold);
    std_threshold = 3;      % from Carrilo-Reid and Jordan Hamm's papers
    static_threshold = .01;
    Spikes = rasterizeDFoF(DeltaFoverF,std_threshold,static_threshold);
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
    spikeTrials = [];
    trialLength = 330;
    figure,
    for i = 1:size(Spikes,2)/trialLength
        spikeTrials{i} = Spikes(:,((i-1)*trialLength+1):i*trialLength);
        Show_Spikes(spikeTrials{i});
%         DeltaTrials(:,:,i) = DeltaFoverF(:,((i-1)*trialLength+1):i*trialLength);
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
    
    

    % Ensemble Analysis
    % figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
    factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
    Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);
    
    % Ensemble stats
    Ensemble = ensembleMetric(Ensemble,AverageImage,ROIcentroid);
    Ensemble = ensembleStat(Ensemble);
    close all
    
    %% Run the analysis again but using the sensory driven response window
    SpikesSen = [];
    for i = 1:length(spikeTrials)
        win = [60 240];
        SpikesSen = horzcat(SpikesSen,spikeTrials{i}(:,win(1):win(2)));
    end
    factorCorrection = 5*floor(size(SpikesSen,2)/5); % Correct for frame size aquisition
    Ensemblesensory = ensembleAnalysis(SpikesSen(:,1:factorCorrection),ROIcentroid);
    
    % Ensemble stats
    Ensemblesensory = ensembleMetric(Ensemblesensory,AverageImage,ROIcentroid);
    Ensemblesensory = ensembleStat(Ensemblesensory);
    close all
    
    %% Save data
    if ~exist([file(fileNum).folder '\output'],'dir')
        mkdir([file(fileNum).folder '\output']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\output'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\output'],[file_name,'.mat']);
    save(savepath,'files', 'Ensemble','spikeTrials','Ensemblesensory','Spikes','ROI', 'ROIcentroid' ,'DeltaFoverF');
    clearvars -except file numFile fileNum filetype foldername
end