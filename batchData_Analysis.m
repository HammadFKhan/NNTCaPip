function batchData = batchData_Analysis
pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
directory = dir(pathname);
L = length(directory);
std_threshold = 2.5;
static_threshold = .7;
f = waitbar(0,'Batch Processing...');
for i = 3:L
     k = i-2;
     waitbar(k/(L-2),f)
     fullfile = strcat(directory(i).folder,'\',directory(i).name);
     load(fullfile);
     Spikes = Spike_Detector_Single(dDeltaFoverF,std_threshold,static_threshold);
     cell_count = length(ROI);
     time = time_adjust(num_images,15);
     corr = correlation_dice(Spikes);
     Connected_ROI = Connectivity_dice(corr, ROI);
     corr_jaccard = correlation_jaccard(Spikes);
    [coactive_cells,detected_spikes] = coactive_index(Spikes);
    calcium_avg = STA(DeltaFoverF,Spikes,std_threshold,10);
    Spikes_shuffled = tempShuffle(Spikes,num_images,cell_count);
    Event_shuffled = spatialShuffle(Spikes,num_images,cell_count);
    surrogate = 5;
    Total_shuffled = allShuffle(Spikes,num_images,cell_count,surrogate);
    shuff_corr = correlation_dice(Total_shuffled);
    [NumActiveNodes,NumNodes,NumEdges,SpatialCentroid,SpatialCentroidVariance,ActivityCentroid,ActivityCentroidVariance]...
        = Network_Analysis(ROIcentroid,Connected_ROI);
    
    [vectorized,sim_index] = cosine_similarity(Spikes,20);
    [shufvectorized,shufsim_index] = cosine_similarity(Total_shuffled,20);
    
    batchData(k).Spikes = Spikes;
    batchData(k).cell_count = length(Spikes(:,1));
    batchData(k).time = time;

    batchData(k).corr = corr;
    batchData(k).Connected_ROI = Connected_ROI;
    batchData(k).corr_jaccard = corr_jaccard;
    batchData(k).coactive_cells = coactive_cells;
    batchData(k).detected_spikes = detected_spikes;
    batchData(k).calcium_avg = calcium_avg;

    batchData(k).Spikes_shuffled = Spikes_shuffled;
    batchData(k).Event_shuffled = Event_shuffled;
    batchData(k).Total_shuffled = Total_shuffled;
    batchData(k).shuff_corr = shuff_corr;
    batchData(k).vectorized = vectorized;
    batchData(k).sim_index = sim_index;
    batchData(k).shufvectorized = shufvectorized;
    batchData(k).shufsim_index = shufsim_index;
    
    
    batchData(k).NumActiveNodes = NumActiveNodes;
    batchData(k).NumNodes = NumNodes;
    batchData(k).NumEdges = NumEdges;
    batchData(k).SpatialCentroid = SpatialCentroid;
    batchData(k).SpatialCentroidVariance = SpatialCentroidVariance;
    batchData(k).ActivityCentroid = ActivityCentroid;
    batchData(k).ActivityCentroidVariAnce = ActivityCentroidVariance;

end
close(f),disp('Batch Process Complete');
