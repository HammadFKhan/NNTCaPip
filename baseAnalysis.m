%% Remove ROIs
if exist('badComponents','var') && ~exist('badComFlag','var')
    [DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A] = ...
        removeROI(DeltaFoverF,dDeltaFoverF,ROI,ROIcentroid,Noise_Power,A,unique(badComponents));
    badComFlag = 1;
end
%% Fix centroids
ROIcentroid = [];
for i = 1:length(ROI)
    if ~isempty(ROI{i})
        blah = vertcat(ROI{i}{:});
        ROIcentroid(i,:) = floor(mean(blah,1));
        continue;
    end
end
%% Analysis
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'));
std_threshold = 2.5;
static_threshold = .01;
Spikes = Spike_Detector_Single(diff(DeltaFoverF),std_threshold,static_threshold);
figure,stack_plot(DeltaFoverF,1,3,1) % Show fluorescence for each cell
figure,Show_Spikes(Spikes) % Plot binary raster plot
%% Ensemble Analysis
% figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);
