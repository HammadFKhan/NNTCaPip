function memMap_CaImAn()
global nam
global memfig
%% setup path to file and package
gcp;                            % start cluster
addpath(genpath('utilities'));
addpath(genpath('deconvolution'));
addpath(genpath('NoRMCorre'));
addpath(genpath('imread_big.m'));

if isfile(nam) == 0
    [filename, pathname] = uigetfile({'*.tiff;*.tif'}, 'Pick a image video file');
    if isequal(filename,0) || isequal(pathname,0)
        disp('User pressed cancel')
    else
        disp(['User selected ', fullfile(pathname, filename)])
    end
    nam = strcat(pathname,filename);          % insert path to tiff stack here
end

sframe=1;

is_memmaped = true;        % choose whether you want to load the file in memory or not
%% Motion correct and memMap Data
%Checks if temporary motion correction file exists (deletes if it does)
% if isfile('motion_corrected.mat') || isfile('motion_corrected.tif')
%     delete('motion_corrected.mat')
%     delete('motion_corrected.tif')
% end
% if exist([nam(1:end-3),'mat'],'file')
%     delete([nam(1:end-3),'mat']);
% end
M2 = memMapmotionCorrection(nam);

%% load file

if exist([nam(1:end-3),'mat'],'file')
    data = matfile([nam(1:end-3),'mat'],'Writable',true);
else
    sframe=1;						% user input: first frame to read (optional, default 1)
    num2read=[];					% user input: how many frames to read   (optional, default until the end)
    chunksize=2000;                 % user input: read and map input in chunks (optional, default read all at once)
    data = memmap_file(nam,sframe,num2read);
    %data = memmap_file_sequence(foldername);
end
% M2 = memMapmotionCorrection(data);

sizY = size(data,'Y');                    % size of data matrix
%% Set parameters
patch_size = [128,128];                   % size of each patch along each dimension (optional, default: [32,32])
overlap = [6,6];                        % amount of overlap in each dimension (optional, default: [4,4])

patches = construct_patches(sizY(1:end-1),patch_size,overlap);
K = 7;                  % number of components to be found per patch
tau = 7;                 % std of gaussian kernel (size of neuron)
p = 2;                   % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
merge_thr = 0.8;         % merging threshold

options = CNMFSetParms(...
    'init_method','greedy',...
    'beta',0.80,...,
    'snmf_max_iter',200,...                     % max # of sparse NMF iterations
    'd1',sizY(1),'d2',sizY(2),...
    'nb',1,...                                  % number of background components per patch
    'gnb',3,...                                 % number of global background components
    'ssub',2,...
    'tsub',4,...
    'p',p,...                                   % order of AR dynamics
    'merge_thr',merge_thr,...                   % merging threshold
    'gSig',tau,...
    'spatial_method','regularized',...
    'cnn_thr',0.2,...
    'patch_space_thresh',0.25,...
    'min_SNR',2,...
    'save_memory',1,...
    'create_memmap',1);

%% Run on patches

[A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data,K,patches,tau,p,options);

%% classify components

rval_space = classify_comp_corr(data,A,C,b,f,options);
ind_corr = rval_space > options.space_thresh;           % components that pass the space correlation test

try  % matlab 2017b or later is needed for the CNN classifier
    [ind_cnn,value] = cnn_classifier(A,[options.d1,options.d2],'cnn_model',options.cnn_thr);
catch
    ind_cnn = true(size(A,2),1);
end

fitness = compute_event_exceptionality(C+YrA,options.N_samples_exc,options.robust_std); % event exceptionality
ind_exc = (fitness < options.min_fitness);

keep = (ind_corr | ind_cnn) & ind_exc;

%% re-estimate temporal components
A_throw = A(:,~keep);
C_throw = C(~keep,:);
A_keep = A(:,keep);
C_keep = C(keep,:);
options.p = 2;      % perform deconvolution
P.p = 2;
[A2,b2,C2] = update_spatial_components(data,C_keep,f,[A_keep,b],P,options);
[C2,f2,P2,S2,YrA2] = update_temporal_components_fast(data,A2,b2,C2,f,P,options);
%% run GUI for modifying component selection (optional, close twice to save values)
Cn = correlation_image_max(data);  % background image for plotting
run_GUI = false;
if run_GUI
    Coor = plot_contours(A,Cn,options,1); close;
    GUIout = ROI_GUI(A,options,Cn,Coor,keep,ROIvars);
    options = GUIout{2};
    keep = GUIout{3};
end
%% plot results
global memfig
memfig = figure(1);
[Coor,json_file] = plot_contours(A2,Cn,options,1);
%     plot_components_GUI(data,A2,C2,b,f2,Cn,options);
%% Output Variables
global AverageImage
global ROIcentroid
global ROI
global num_images
global DeltaFoverF
global dDeltaFoverF
global Noise_Power
global files
global A
global C
global ops
A = A2;
C = Cn;
ops = options;
files = nam;
[H,W] = size(Cn);
AverageImage =ones(H,W);
num_images = size(C2,1);
DeltaFoverF = C2;
dDeltaFoverF = S2;
ROIcentroid = {json_file.centroid}';ROIcentroid = cell2mat(ROIcentroid);
ROI = {};
for i = 1:length(Coor)
    for ii = 1:length(Coor{i})
        ROI{i,1}{ii,1} = Coor{i}(1:2,ii)';
    end
end
Noise_Power = P;
end