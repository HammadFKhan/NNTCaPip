%% Run the base pipeline for denoised/deconvolved calcium traces
%% CaimAn Single File ROI Extraction (Large Dataset)
clear
clc
close all;
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'));
addpath(genpath('Pipelines'))
global nam
global memfig
batch = 1;
if batch == 1
    pathname = strcat(uigetdir(pwd,'Input Directory'),'\');
    savepathname = strcat(uigetdir(pwd,'Output Directory'),'\');
    directory = dir(fullfile(pathname,'*.tif'));
    L = length(directory);
    for idx = 1:L
        clear files AverageImage num_images...
            DeltaFoverF dDeltaFoverF ROIcentroid ROI Noise_Power
        filename = directory(idx).name
        nam = strcat(pathname,filename);
        tic
        Start_MemMap_CaImAn
        toc
        savepath = strcat(savepathname,filename,'.mat');
        save(savepath,'files','AverageImage','num_images',...
            'DeltaFoverF','dDeltaFoverF','ROIcentroid','ROI','Noise_Power','C','A','ops');
        try
            savepathfig = strcat(savepathname,filename(1:end-4),'.fig');
            saveas(memfig,savepathfig);
        catch ME
            warning('Contour figure not saved')
            continue
        end
        disp('Saved!')
    end
else
    nam = '';
    tic
    Start_MemMap_CaImAn
    toc
end
