%% CaimAn File ROI Extraction for large datasets
% Run the base pipeline for denoised/deconvolved calcium traces
% Folder directory is chosen by the user
% All Tiff files in the subdirectory are motion corrected and stored as
% .h5 scientific files for memory mapping large data sets. 
% NOTE: motion correction is skipped if it has already been done
% All motion corrected data is denoised and deconvolved using OASIS.

% Pre_shifts_nr.mat provides motion corrected data
% ds_data provides memmap data for pointing certain parts of the pipeline.
% NOTE: ds_data allocation is skipped if it already has been done
% All variables for analysis are stored under output

clear
clc
close all;
set(0,'DefaultFigureWindowStyle','normal')
addpath(genpath('main'))
addpath(genpath('Pipelines'))
batchFlag = 0; % Sets if data is to be processed in batches (almost always)
CaImAnFull(batchFlag); 