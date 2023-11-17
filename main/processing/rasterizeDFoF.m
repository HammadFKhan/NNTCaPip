function [Spikes] = rasterizeDFoF(DeltaFoverF,threshold,spikemin)

%% Description 
% % Inputs
% 1. DeltaFoverF - deconvolved spike rate corresponding (s from the CaImAn) dF/F traces (cells x time)
% 2. threshold - standard deviation about the mean of dF/F traces, above which spike is detected 
% 3. spikemin - the bare minimum dF/F value for spike detection threshold (help in rejecting bad/mislabelled cells)

% % Output
% Spikes - 2D matrix of cell x time with 1 representing spike detected at
% the time step

%% Reference - https://github.com/slceto/CalFDR/blob/master/CalFDR.m
%% Reference - DOI: 10.1126/science.aaf7560 - Imprinting and recalling cortical ensembles, Science 2016 - Luis Carrillo-Reid , R. Yuste
%% Reference - https://doi.org/10.1523/JNEUROSCI.5214-14.2015; Carrillo-Ried - Endogenous Sequential Cortical Activity Evoked by Visual Stimuli
Spikes = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2));
derDeltaFoverF = zeros(size(DeltaFoverF,1),size(DeltaFoverF,2));
derDeltaFoverF(:,2:end) = diff(DeltaFoverF,1,2);

noisePerCell = std(derDeltaFoverF,0,2);
threshPerCell = threshold*noisePerCell; % Genreally 2.5 to 3 is choosen according to Luis Carrillo-Ried papers
threshPerCell(threshPerCell < spikemin) = spikemin;
Spikes = derDeltaFoverF - threshPerCell;
Spikes(Spikes>0) = 1;
Spikes(Spikes<=0) = 0; 

end

