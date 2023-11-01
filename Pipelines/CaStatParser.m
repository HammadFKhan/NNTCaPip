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
totalActivity = [];
for i = 1:5
    SDbuffer = correlation_dice(S(i).FOV.Ensemble.ensemble);
    Simbuffer = S(i).FOV.Ensemble.sim_index;
    rCorr(i,:) = [mean(SDbuffer,'all'), std(tril(SDbuffer,-1),[],'all')];
    %     figure,imagesc(SDbuffer),colormap(hot),caxis([0 max(tril(SDbuffer,-1),[],'all')/1.5]),colorbar
    %     figure,imagesc(Simbuffer),colormap(hot),caxis([0 max(tril(Simbuffer,-1),[],'all')/1.5]),colorbar
    %     figure,histogram(SDbuffer(SDbuffer>0.1),'DisplayStyle','stairs'),title(['N: ' num2str(i)])
    blah = tril(SDbuffer,-1);
    SDactivity = [SDactivity;blah(blah>0)];
    Simactivity = [Simactivity;mean(Simbuffer)'];
    
    mean(SDbuffer(SDbuffer>0))
    
    sparseSimactivity = [sparseSimactivity;mean(Simbuffer(Simbuffer>0))];
    blah = tril(SDbuffer,-1);
    sparseSDactivity = [sparseSDactivity;mean(SDbuffer(SDbuffer>0))]; %include nonzero matrix analysis as well
    totalActivity(i) = size(SDbuffer,1);
end
%%
figure,customBoxplot(SDactivity),title('SD activity'),ylim([0 .4]),box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
figure,customBoxplot(sparseSDactivity),title('Sparse SD activity'),ylim([0 .4]);box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
figure,customBoxplot(Simactivity),title('Sim activity'),ylim([0 .2]);box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
figure,customBoxplot(sparseSimactivity),title('Sparse Sim activity'),ylim([0 .4]);box off,set(gca,'TickDir','out'),set(gca,'FontSize',16)
%% Compute ensemble stability
totalEnsembleRecruitment = []; totalMaxRecruitment = []; totalNumEdges = []; totalEnsembleSize = []; totalEnsembleNum = []; totalNumNodes = []; totalEnsembleWeight = [];
for i = 6:9
    Ensemble = ensembleStat(S(i).FOV.Ensemblesensory);
    totalEnsembleRecruitment = [totalEnsembleRecruitment; Ensemble.ensembleRecruitment];
    totalMaxRecruitment = [totalMaxRecruitment; Ensemble.maxRecruitment];
    totalNumEdges = [totalNumEdges;cell2mat(Ensemble.NumEdges)'];
    totalEnsembleSize = [totalEnsembleSize; Ensemble.ensembleSize'];
    totalEnsembleNum = [totalEnsembleNum; floor(Ensemble.ensembleNum)/4];
    totalNumNodes = [totalNumNodes; max(cell2mat(Ensemble.NumNodes))];
    com = vertcat(Ensemble.Connected_ROI{:});
    totalEnsembleWeight = [totalEnsembleWeight;com(:,3)];
end
%%
figure,customBoxplot(totalEnsembleRecruitment),title('Ensemble Recruitment');
figure,customBoxplot(totalMaxRecruitment),title('Max Recruitment');box off, ylim([0.0 .4])
figure,customBoxplot([totalEnsembleRecruitment;totalMaxRecruitment]),title('Both');box off, ylim([0.0 .4])
figure,customBoxplot(totalEnsembleNum),box off,set(gca,'TickDir','out'),title('Ensemble Number')
figure,customBoxplot(totalEnsembleSize),box off,set(gca,'TickDir','out'),title('Ensemble Size')
figure,customBoxplot(totalNumEdges),box off,set(gca,'TickDir','out'),title('Connections')
figure,customBoxplot(totalNumNodes),box off,set(gca,'TickDir','out'),title('Nodes')
figure,customBoxplot(totalEnsembleWeight),box off,set(gca,'TickDir','out'),title('Ensemble Weight')

%% Compute entropy stats