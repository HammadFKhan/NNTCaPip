# NNTCaPip
Image and data processing pipeline for calcium imaging

The main focus of this pipline is to create a one stop shop for analysis of in vivo calcium images. This data pipline includes general fluorescence mapping, spike sourting, and population spike rates. Along with more intensive coactive indexing, spatial temporal correlations, shuffling, and principle component analysis using singular value decomposition for network assemblies. More elaborate explainations are found later in this document. 

## Installation
Please download the full repository into your working MATLAB folder. You must also download the reposity of [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB) for ROI extraction and motion correction of the image stack. You will need to add the [CaImAn](https://github.com/flatironinstitute/CaImAn-MATLAB) folder to the NNTCaPip MATLAB path.

## Running your first calcium video
To run the program open the matlab (.m) file labeled _CONSOLE_. This file acts as a central engine for the entire pipeline (i.e all function calls are eventually mapped to this file). I recommend running the file section by section as the current build has both the single and batch processing in the same space -- to be fixed later. 

### Single vs Batch Analysis
The current build of this pipeline comes with the options of single file or batch processing. 




 _Neurons.tiff_ is an acceptable image stack example provided to test the pipeline. 
