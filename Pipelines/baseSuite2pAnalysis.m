
%% Fix centroids
ROIcentroid = cellfun(@(x) x.med, stat,'UniformOutput',false);
ROIcentroid = vertcat(ROIcentroid{:});

%%% Clean up soma
skew = cellfun(@(x) x.skew, stat);
keepId = skew>1;
Fnew = F(keepId,:);
Fneun = Fneu(keepId,:);
ROIcentroid = ROIcentroid(keepId,:);
%%% Run convolved traces after loading in suite2P data

dF = convolveFluroescence(Fnew,Fneun);
%% Analysis
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'));
std_threshold = 2.5;
static_threshold = .01;
Spikes = rasterizeDFoF(diff(dF),std_threshold,static_threshold);
%figure,stack_plot(dF,1,3,1) % Show fluorescence for each cell
figure,Show_Spikes(Spikes) % Plot binary raster plot
%% Ensemble Analysis
% figure,[Coor,json_file] = plot_contours(A,C,ops,0); % contour plot of spatial footprints
factorCorrection = 5*floor(size(Spikes,2)/5); % Correct for frame size aquisition
Ensemble = ensembleAnalysis(Spikes(:,1:factorCorrection),ROIcentroid);
