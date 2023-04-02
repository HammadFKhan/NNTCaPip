%% Load in calcium data
% clear
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
Simactivity = []; sparseSimactivity = [];

for i = 1:L
    SDbuffer = correlation_dice(S(i).FOV.Ensemble.ensemble);
    Simbuffer = S(i).FOV.Ensemble.sim_index;
    rCorr(i,:) = [mean(SDbuffer,'all'), std(tril(SDbuffer,-1),[],'all')];
%     figure,imagesc(SDbuffer),colormap(hot),caxis([0 max(tril(SDbuffer,-1),[],'all')/1.5]),colorbar
%     figure,imagesc(Simbuffer),colormap(hot),caxis([0 max(tril(Simbuffer,-1),[],'all')/1.5]),colorbar
%     figure,histogram(SDbuffer(SDbuffer>0.1),'DisplayStyle','stairs'),title(['N: ' num2str(i)])
    SDactivity = [SDactivity;mean(SDbuffer)'];
    Simactivity = [Simactivity;mean(Simbuffer)'];
    sparseSimactivity = [Simactivity;mean(Simbuffer(Simbuffer>0))];
    sparseSDactivity = [sparseSDactivity;mean(SDbuffer(SDbuffer>0))]; %include nonzero matrix analysis as well
end
figure,customBoxplot(SDactivity),title('SD activity')
figure,customBoxplot(sparseSDactivity),title('Sparse SD activity');
figure,customBoxplot(Simactivity),title('Sim activity')
figure,customBoxplot(sparseSimactivity),title('Sparse Sim activity')
%% Compute ensemble stability
totalEnsembleRecruitment = []; totalMaxRecruitment = []; totalNumEdges = []; totalEnsembleSize = []; totalEnsembleNum = []; totalNumNodes = [];
for i = 1:L
    Ensemble = ensembleStat(S(i).FOV.Ensemble);
    totalEnsembleRecruitment = [totalEnsembleRecruitment; Ensemble.ensembleRecruitment];
    totalMaxRecruitment = [totalMaxRecruitment; Ensemble.maxRecruitment];
    totalNumEdges = [totalNumEdges;cell2mat(Ensemble.NumEdges)'];
    totalEnsembleSize = [totalEnsembleSize; Ensemble.ensembleSize'];
    totalEnsembleNum = [totalEnsembleNum; Ensemble.ensembleNum];
    totalNumNodes = [totalNumNodes; max(cell2mat(Ensemble.NumNodes))];
end
%%
figure,customBoxplot(totalEnsembleRecruitment),title('Ensemble Recruitment');
figure,customBoxplot(totalMaxRecruitment),title('Max Recruitment');box off, ylim([0.0 .4])
figure,customBoxplot([totalEnsembleRecruitment;totalMaxRecruitment]),title('Both');box off, ylim([0.0 .4])
figure,customBoxplot(totalEnsembleNum/4),box off,set(gca,'TickDir','out'),title('Ensemble Number')
figure,customBoxplot(totalEnsembleSize),box off,set(gca,'TickDir','out'),title('Ensemble Size')
figure,customBoxplot(totalNumEdges),box off,set(gca,'TickDir','out'),title('Connections')
figure,customBoxplot(totalNumNodes),box off,set(gca,'TickDir','out'),title('Nodes')

%% Compute entropy stats