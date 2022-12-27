%% Load in calcium data
clear
pathname = uigetdir(pwd,'Input Directory');
pathname = fullfile(pathname);
directory = dir(fullfile(pathname,'*.mat'));
L = length(directory);
for i = 1:L
    disp(['Parsing: ' num2str(directory(i).name)])
    S(i).FOV = load(fullfile(directory(i).folder,directory(i).name));
end 
%% Compute and concatenate sorenson dice pairwise correlation matrix across animals/FOVs
SDactivity = [];
sparseSDactivity = [];
for i = 1:L
    SDbuffer = correlation_dice(S(i).FOV.Ensemble.ensemble);
    SDactivity = [SDactivity;mean(SDbuffer)'];
    sparseSDactivity = [sparseSDactivity;mean(SDbuffer(SDbuffer>0))]; %include nonzero matrix analysis as well
end
figure,boxplot(SDactivity),title('SD activity')
figure,boxplot(sparseSDactivity),title('Sparse SD activity');ylim([0.04 0.3])
%% Compute ensemble stability
totalEnsembleRecruitment = []; totalMaxRecruitment = [];
for i = 1:L
    totalEnsembleRecruitment = [totalEnsembleRecruitment; S(i).FOV.Ensemble.ensembleRecruitment];
    totalMaxRecruitment = [totalMaxRecruitment; S(i).FOV.Ensemble.maxRecruitment];
end
figure,boxplot(totalEnsembleRecruitment),title('Ensemble Recruitment');
figure,boxplot(totalMaxRecruitment),title('Max Recruitment');box off, ylim([0.0 .4])
figure,boxplot([totalEnsembleRecruitment;totalMaxRecruitment]),title('Both');box off, ylim([0.0 .4])
%% Compute entropy stats