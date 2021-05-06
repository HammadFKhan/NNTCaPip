# NNTCaPip

![Neurons](https://user-images.githubusercontent.com/56140216/117201664-ea381b80-adba-11eb-9ffe-b2bed4401cf5.png)

### Image and data processing pipeline for calcium imaging

The main focus of this pipline is to create a one stop shop for analysis of in vivo calcium images. This data pipline includes general fluorescence mapping, spike sorting, and population spike rates. Along with more intensive coactive indexing, spatial temporal correlations, shuffling, and principle component analysis using singular value decomposition for network assemblies. More elaborate explainations are found later in this document. 

## Installation

### Requirements

- MATLAB 2016     
- Parallel Computing Toolbox
- Signal Processing Toolbox

Please download the full repository into your working MATLAB folder. You must also download the reposity of [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB) for ROI extraction and motion correction of the image stack. You will need to add the [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB) folder to the NNTCaPip MATLAB path.

## Running your first calcium video
To run the program open the matlab (.m) file labeled _CONSOLE_. This file acts as a central engine for the entire pipeline (i.e all function calls are eventually mapped to this file). I recommend running the file section by section as the current build has both the single and batch processing in the same space -- to be fixed later. 

### Single vs Batch Analysis
The current build of this pipeline comes with the options of single file or batch processing. To run either a single tiff stack or example video provided in the pipeline run the MATLAB section labeled _Caiman Single File Analysis_ followed by _Single File Analysis_. Both these section process the bulk of the data. Now you are ready to visualize all of it! To do this, run the section labeled _Plot all_. 

```matlab
%% CaimAn Single Batch Analysis
clear
clc
close all;
tic
Start_CaImAn
toc
```
```matlab
%% Single File Analysis
set(0,'DefaultFigureWindowStyle','docked')...
```
```matlab
%% Plot all the Figures
addpath('Figures');...
```

## Understanding the Data

All fluorescence data is analyzed using a ROI x Frame matrix. Further, the bulk of spike sorting analysis is done by converting arbitrary fluorescence values into a binary matrix. To understand what the data is telling you please refer to our [Wiki](https://github.com/Neurohm/NNTCaPip/wiki) page.
