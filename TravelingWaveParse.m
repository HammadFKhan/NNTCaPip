pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S(i).FOV = load(fullfile(directory(i).folder,directory(i).name));
end 
%%
totalActivityCentroid = horzcat(totalActivityCentroid,S(4).FOV.ActivityCentroid);
totalConnectedROI = horzcat(totalConnectedROI,S(4).FOV.Connected_ROI);
totalROIcentroid = vertcat(totalROIcentroid, S(4).FOV.ROIcentroid);
totalSpikes = horzcat(S(1).FOV.Spikes,S(2).FOV.Spikes);
%%
figure,
NodeSize = 1;EdgeSize = 1;
for i = 1:length(lateSpikeTrials(lateSpikeTrials<123))
    subplot(9,9,i),Cell_Map_Dice(AverageImage,totalConnectedROI{lateSpikeTrials(i)},totalROIcentroid,NodeSize,EdgeSize)
end
suptitle('Late Spike')

figure,
NodeSize = 1;EdgeSize = 1;
for i = 1:length(nolateSpikeTrials(nolateSpikeTrials<123))
    subplot(6,6,i),Cell_Map_Dice(AverageImage,totalConnectedROI{nolateSpikeTrials(i)},totalROIcentroid,NodeSize,EdgeSize)
end
suptitle('No Late Spike')
