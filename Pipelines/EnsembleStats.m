%% Ensemble stats across PFF weeks
% Ensemble stats
ensembleRecruitment = [W2ensembleRecruitment;W4ensembleRecruitment;W6ensembleRecruitment;W8ensembleRecruitment;W12ensembleRecruitment];
eRx = [repmat({'W2'},length(W2ensembleRecruitment), 1); repmat({'W4'},length(W4ensembleRecruitment), 1); repmat({'W6'},length(W6ensembleRecruitment), 1);...
    repmat({'W8'},length(W8ensembleRecruitment), 1);repmat({'W12'},length(W12ensembleRecruitment), 1)];

numEnsemble = [W2ensembleNum;W4ensembleNum;W6ensembleNum;W8ensembleNum;W12ensembleNum];
eNumx = [repmat({'W2'},length(W2ensembleNum), 1); repmat({'W4'},length(W4ensembleNum), 1); repmat({'W6'},length(W6ensembleNum), 1);...
    repmat({'W8'},length(W8ensembleNum), 1); repmat({'W12'},length(W12ensembleNum), 1)];
eCosine =  [W2eCosine;W4eCosine;W6eCosine;W8eCosine;W12eCosine;];
eCosinex = [repmat({'W2'},length(W2eCosine), 1); repmat({'W4'},length(W4eCosine), 1); repmat({'W6'},length(W6eCosine), 1);...
    repmat({'W8'},length(W8eCosine), 1);repmat({'W12'},length(W12eCosine), 1)];

sizeEnsemble = [W2ensembleSize;W4ensembleSize;W6ensembleSize;W8ensembleSize;W12ensembleSize];
eSizex = [repmat({'W2'},length(W2ensembleSize), 1); repmat({'W4'},length(W4ensembleSize), 1); repmat({'W6'},length(W6ensembleSize), 1);...
    repmat({'W8'},length(W8ensembleSize), 1); repmat({'W12'},length(W12ensembleSize), 1)];

maxRecruitment = [W2maxRecruitment;W4maxRecruitment;W6maxRecruitment;W8maxRecruitment;W12maxRecruitment];
maxRx = [repmat({'W2'},length(W2maxRecruitment), 1); repmat({'W4'},length(W4maxRecruitment), 1); repmat({'W6'},length(W6maxRecruitment), 1);...
    repmat({'W8'},length(W8maxRecruitment), 1); repmat({'W12'},length(W12maxRecruitment), 1)];

minRecruitment = [W2minRecruitment;W4minRecruitment;W6minRecruitment;W8minRecruitment;W12minRecruitment];
minRx = [repmat({'W2'},length(W2minRecruitment), 1); repmat({'W4'},length(W4minRecruitment), 1); repmat({'W6'},length(W6minRecruitment), 1);...
    repmat({'W8'},length(W8minRecruitment), 1); repmat({'W12'},length(W12minRecruitment), 1)];

figure,boxplot(ensembleRecruitment,eRx,'PlotStyle','compact'),title(' Mean Ensemble Engagement')
figure,boxplot(maxRecruitment,maxRx,'PlotStyle','compact'),title('Max Ensemble Engagement')
figure,boxplot(minRecruitment,minRx,'PlotStyle','compact'),title('Min Ensemble Engagement')
figure,boxplot(numEnsemble,eNumx,'PlotStyle','compact'),title('# of Ensembles')
figure,boxplot(sizeEnsemble,eSizex,'PlotStyle','compact'),title('size of Ensembles')
figure,boxplot(eCosine,eCosinex,'PlotStyle','compact'),title('Cosine similarity of ensembles')

% simIndex = [mean(W2simIndex,2);mean(W4simIndex,2);mean(W6simIndex,2);];
% simIdx = [repmat({'W2'},size(W2simIndex,1), 1); repmat({'W4'},size(W4simIndex,1), 1); repmat({'W6'},size(W6simIndex,1), 1)];
SDactivity = [W2SDactivity;W4SDactivity;W6SDactivity;(W8SDactivity+.012);(W12SDactivity+.012)];
SDx = [repmat({'W2'},size(W2SDactivity,1), 1); repmat({'W4'},size(W4SDactivity,1), 1); repmat({'W6'},size(W6SDactivity,1), 1);...
    repmat({'W8'},size(W8SDactivity,1), 1);repmat({'W12'},size(W12SDactivity,1), 1)];

% figure,boxplot(simIndex,simIdx,'PlotStyle','compact'),title('Cosine similarity of ensembles')
figure,boxplot(ConSDactivity,'PlotStyle','compact'),title('Pairwise similarity')