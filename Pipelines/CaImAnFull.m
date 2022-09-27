function CaImAnFull(batchFlag)
%% setup path to file and package
gcp;                            % start cluster
addpath(genpath('CaImAn-MATLAB'))
addpath(genpath('utilities'));
addpath(genpath('deconvolution'));
addpath(genpath('NoRMCorre')); %add the NoRMCorre motion correction package to MATLAB path

% foldername = 'E:\In vivo Imaging\09242022-Rbp4_23773\Calcium\preprocessed\Test\';
%          % folder where all the files are located.
% filetype = 'tif'; % type of files to be processed
%         % Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
% file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)

if batchFlag
    foldername = strcat(uigetdir(pwd,'Input Directory'),'\');
    filetype = 'tif'; % type of files to be processed
    % Types currently supported .tif/.tiff, .h5/.hdf5, .raw, .avi, and .mat files
    file = subdir(fullfile(foldername,['*.',filetype]));   % list of filenames (will search all subdirectories)
else
    [filename, foldername] = uigetfile({'*.tiff;*.tif'}, 'Pick a image video file');
    if isequal(filename,0) || isequal(foldername,0)
        disp('User pressed cancel')
    else
        disp(['User selected ', fullfile(foldername, filename)])
    end
    file = subdir(fullfile(foldername,filename));
end


FOV = size(read_file(file(1).name,1,1));
numFiles = length(file);

%% motion correct (and save registered h5 files as 2d matrices (to be used in the end)..)
% register files one by one. use template obtained from file n to
% initialize template of file n + 1;

motion_correct = true;                            % perform motion correction
non_rigid = true;                                 % flag for non-rigid motion correction
output_type = 'h5';                               % format to save registered files

if non_rigid; append = '_nr'; else; append = '_rig'; end        % use this to save motion corrected files

for fileNum = 1:numFiles
    options_mc = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[128,128],'init_batch',200,...
        'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
        'use_parallel',true,'output_type',output_type);
    
    template = [];
    col_shift = [];
    fullname = file(fileNum).name;
    [folder_name,file_name,ext] = fileparts(fullname);
    output_filename = fullfile(folder_name,[file_name,append,'.',output_type]);
    options_mc = NoRMCorreSetParms(options_mc,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name
    if ~exist(output_filename,'file') %check is motion correct file already exists
        if motion_correct
            [M,shifts,template,options_mc,col_shift] = normcorre_batch(fullname,options_mc,template);
            save(fullfile(folder_name,[file_name,'_shifts',append,'.mat']),'shifts','-v7.3');           % save shifts of each file at the respective folder
        else    % if files are already motion corrected convert them to h5
            convert_file(fullname,'h5',fullfile(folder_name,[file_name,'_mc.h5']));
        end
    end
    
end

if motion_correct
    %     registered_files = subdir(fullfile(foldername,['*',append,'.',output_type]));  % list of registered files (modify this to list all the motion corrected files you need to process)
    registered_files = subdir(fullfile(foldername,['*',append,'.',output_type]));
    if fileNum==1 % Checks to see if we only wanted to look at one file
        registered_files = subdir(fullfile(foldername,[file_name,append,'.',output_type]));
    end
else
    registered_files = subdir(fullfile(foldername,'*_mc.h5'));
end
%% Big loop for all files
for fileNum = 1:numFiles
    fr = 30;                                         % frame rate
    tsub = 8;                                        % degree of downsampling (for 30Hz imaging rate you can try also larger, e.g. 8-10)
    ds_filename = [foldername,'/',file_name,'_ds_data.mat'];
    if ~exist(ds_filename,'file') %check if downsampled data file already exists
        data = matfile(ds_filename,'Writable',true);
        [~,fileInfo] = bigread2(file(fileNum).name); % Read in the tif file just for indexing
        data.fileInfo = fileInfo;
        data_type = fileInfo.dataType;
        FOV = fileInfo.FOV;
        dims = FOV;
        ndimsY = length(dims);                       % number of dimensions (data array might be already reshaped)
        Ts = dims(end);
        F_dark = Inf;                                    % dark fluorescence (min of all data)
        batch_size = 2000;                               % read chunks of that size
        batch_size = round(batch_size/tsub)*tsub;        % make sure batch_size is divisble by tsub
        % Ts = zeros(numFiles,1);                        % store length of each file
        cnt = 0;                                         % number of frames processed so far
        tt1 = tic;
        
        name = registered_files(fileNum).name;
        info = h5info(name);
        Ysub = zeros(FOV(1),FOV(2),floor(Ts/tsub),data_type);
        cnt_sub = 0;
        for t = 1:batch_size:Ts
            Y = read_file(name,t,min(batch_size,Ts-t+1));
            F_dark = min(nanmin(Y(:)),F_dark);
            ln = size(Y,ndimsY);
            %     Y = reshape(Y,[FOV,ln]);
            Y = cast(downsample_data(Y,'time',tsub),data_type);
            ln = size(Y,3);
            Ysub(:,:,cnt_sub+1:cnt_sub+ln) = Y;
            cnt_sub = cnt_sub + ln;
        end
        data.Y= Ysub;
        data.Yr = reshape(Ysub,[],cnt_sub);
        toc(tt1);
        % cnt = cnt + cnt_sub;
        data.sizY = size(Ysub);
        data.F_dark = F_dark;
    else
        data = matfile(ds_filename,'Writable',true);
        try
            fileInfo = data.fileInfo;
        catch % Future build correction if fileInfo did not exist in data
            [~,fileInfo] = bigread2(file(fileNum).name);
            data.fileInfo = fileInfo;
        end
        FOV = fileInfo.FOV;
        dims = FOV;
        Ts = dims(end);
    end
    
    
    
    %% Set parameters
    sizY = data.sizY;                       % size of data matrix
    patch_size = [128,128];                   % size of each patch along each dimension (optional, default: [128,128])
    overlap = [6,6];                        % amount of overlap in each dimension (optional, default: [6,6])
    
    patches = construct_patches(sizY(1:end-1),patch_size,overlap);
    K = 25;                  % number of components to be found per patch
    tau = 4;                 % std of gaussian kernel (size of neuron)
    p = 2;                   % order of autoregressive system (p = 0 no dynamics, p=1 just decay, p = 2, both rise and decay)
    merge_thr = 0.4;         % merging threshold
    
    options = CNMFSetParms(...
        'init_method','sparse_NMF',...              % Segmentation type ('greedy' for soma 'sparse_NMF' for dendrites
        'beta',0.5,...,                             % NMF converging coefficient (higher is stricter) (default:0.9)
        'snmf_max_iter',80,...                      % max # of sparse NMF iterations (default:50)
        'd1',sizY(1),'d2',sizY(2),...               % FOV size (512x512 typically)
        'nb',1,...                                  % number of background components per patch
        'gnb',2,...                                 % number of global background components
        'ssub',1,...                                % spatial downsample (will reupsample) (default:1)
        'tsub',2,...                                % temporal downsample (will reupsample) (default:2)
        'fr',fr/tsub,...                            % downsamples
        'p',p,...                                   % order of AR dynamics
        'merge_thr',merge_thr,...                   % merging threshold (default:0.4)
        'gSig',tau,...                              % body size (default:0.4)
        'spatial_method','regularized',...          % spatial threshold
        'cnn_thr',0.2,...                           % classifier threshold (default:0.2)
        'patch_space_thresh',0.25,...               % merge patch threshold
        'min_SNR',2,...                             % minimum signal SNR
        'search_method','dilate');                 % method for determining footprint of spatial components 'ellipse' or 'dilate' (default: 'dilate')
    
    %% Run on patches
    
    [A,b,C,f,S,P,RESULTS,YrA] = run_CNMF_patches(data.Y,K,patches,tau,0,options);  % do not perform deconvolution here since we have downsampled
    
    %% compute correlation image on a small sample of the data (optional - for visualization purposes)
    Cn = correlation_image_max(data,8);
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
    
    %% view contour plots of selected and rejected components (optional)
    throw = ~keep;
    Coor_k = [];
    Coor_t = [];
    
    
    %% keep only the active components
    
    A_keep = A(:,keep);
    A_throw = A(:,throw);
    C_keep = C(keep,:);
    
    %% extract residual signals for each trace
    
    if exist('YrA','var')
        R_keep = YrA(keep,:);
    else
        R_keep = compute_residuals(data,A_keep,b,C_keep,f);
    end
    %% extract fluorescence on native temporal resolution
    
    options.fr = options.fr*tsub;                   % revert to origingal frame rate
    N = size(C_keep,1);                             % total number of components
    T = sum(Ts);                                    % total number of timesteps
    C_full = imresize(C_keep,[N,T]);                % upsample to original frame rate
    R_full = imresize(R_keep,[N,T]);                % upsample to original frame rate
    F_full = C_full + R_full;                       % full fluorescence
    f_full = imresize(f,[size(f,1),T]);             % upsample temporal background
    
    S_full = zeros(N,T);
    
    chunkSize = 2000;
    P.p = 0;
    ind_T = [0;cumsum(Ts(:))];
    options.nb = options.gnb;
    % for i = 1:numFiles
    %     inds = ind_T(i)+1:ind_T(i+1);   % indeces of file i to be updated
    %     [C_full(:,inds),f_full(:,inds),~,~,R_full(:,inds)] = update_temporal_components_fast(registered_files(i).name,A_keep,b,C_full(:,inds),f_full(:,inds),P,options);
    %     disp(['Extracting raw fluorescence at native frame rate. File ',num2str(i),' out of ',num2str(numFiles),' finished processing.'])
    % end
    
    %% extract DF/F and deconvolve DF/F traces
    
    [F_dff,F0] = detrend_df_f(A_keep,b,C_full,f_full,R_full,options);
    
    C_dec = zeros(N,T);         % deconvolved DF/F traces
    S_dec = zeros(N,T);         % deconvolved neural activity
    bl = zeros(N,1);            % baseline for each trace (should be close to zero since traces are DF/F)
    neuron_sn = zeros(N,1);     % noise level at each trace
    g = cell(N,1);              % discrete time constants for each trace
    if p == 1; model_ar = 'ar1'; elseif p == 2; model_ar = 'ar2'; else; error('This order of dynamics is not supported'); end
    
    for i = 1:N
        spkmin = options.spk_SNR*GetSn(F_dff(i,:));
        lam = choose_lambda(exp(-1/(options.fr*options.decay_time)),GetSn(F_dff(i,:)),options.lam_pr);
        [cc,spk,opts_oasis] = deconvolveCa(F_dff(i,:),model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
            'window',150,'lambda',lam,'smin',spkmin);
        bl(i) = opts_oasis.b;
        C_dec(i,:) = cc(:)' + bl(i);
        S_dec(i,:) = spk(:);
        neuron_sn(i) = opts_oasis.sn;
        g{i} = opts_oasis.pars(:)';
        disp(['Performing deconvolution. Trace ',num2str(i),' out of ',num2str(N),' finished processing.'])
    end
    %% Plot
    % Save Data for video plotting
    data.b = b;data.f = f;data.A = A_keep; data.C = C_keep;data.R = R_keep;
    savfig = figure(1);
    ax1 = subplot(121); [Coor,json_file] = plot_contours(A_keep,Cn,options,0); title('Selected components','fontweight','bold','fontsize',14);
    ax2 = subplot(122); plot_contours(A_throw,Cn,options,0);title('Rejected components','fontweight','bold','fontsize',14);
    linkaxes([ax1,ax2],'xy')
    %% Output Variables
    % global AverageImage
    % global ROIcentroid
    % global ROI
    % global num_images
    % global DeltaFoverF
    % global dDeltaFoverF
    % global Noise_Power
    % global files
    % global A
    % global C
    % global ops
    A = A_keep; %spatial components
    C = Cn;
    ops = options;
    files = file(fileNum).name;
    [H,W] = size(Cn);
    AverageImage =ones(H,W);
    num_images = FOV(3);
    DeltaFoverF = C_dec; %deconvolved calcium traces
    dDeltaFoverF = S_dec; % deconvolved spike trains
    ROIcentroid = {json_file.centroid}';ROIcentroid = cell2mat(ROIcentroid); % ROI maps for ensembles
    ROI = {};
    for i = 1:length(Coor)
        for ii = 1:length(Coor{i})
            ROI{i,1}{ii,1} = Coor{i}(1:2,ii)';
        end
    end
    Noise_Power = R_full;
    
    if ~exist([file(fileNum).folder '\output'],'dir')
        mkdir([file(fileNum).folder '\output']);
    end
    [folder_name,file_name,~] = fileparts(file(fileNum).name);
    if exist(fullfile([folder_name, '\output'],[file_name,'.mat']),'file')
        file_name = [file_name '_' datestr(now,30) '_'];
    end
    savepath = fullfile([folder_name, '\output'],[file_name,'.mat']);
    save(savepath,'files','AverageImage','num_images',...
        'DeltaFoverF','dDeltaFoverF','ROIcentroid','ROI','Noise_Power','C','A','ops');
    try
        savepathfig = fullfile([folder_name, '\output'],[file_name,'.fig']);
        saveas(savfig,savepathfig);
        close(savfig)
    catch ME
        warning('Contour figure not saved')
    end
    disp('Saved!')
end
end